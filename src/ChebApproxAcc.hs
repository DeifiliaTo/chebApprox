{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DeriveAnyClass  #-}
{-# LANGUAGE DeriveGeneric   #-}
{-# LANGUAGE PatternSynonyms #-}

module  ChebApproxAcc
where
    import Prelude as P
    import Data.Array.Accelerate as A
    import Data.Array.Accelerate.Debug as A
    import Data.Array.Accelerate.Interpreter as I
    import Data.Array.Accelerate.LLVM.Native as CPU


   
    value :: Exp (Int)
    value = (constant 2)

    value2 :: Exp (Int)
    value2 = (constant 1)   

    value3 :: Exp (Double)
    value3 = (constant 0.02)

    arr1 :: Acc (Vector Double)
    arr1 = use (fromList (Z :. 10) [0..])

    vec0 :: Acc (Vector Double)
    vec0 = use (fromList (Z :. 1) [1..])

    vec1 :: Acc (Vector Double)
    vec1 = use (fromList (Z :. 2) [0..])

    -- Returns value of chebyshev zero
    computeChebNode :: Exp (Int) -> Exp (Int) -> Exp (Double)
    computeChebNode n k = 
        cos ((2*A.fromIntegral(n)+1-2*A.fromIntegral(k))*pi /(2*A.fromIntegral(n)+2))
  
     -- Creates list of chebyshev zeros
    chebNodes :: Exp Int -> Acc (Vector Double)
    chebNodes n =
         let nodes = enumFromN (lift (Z:.(n+1))) 0 in
            A.map (\x -> computeChebNode (n) x) (nodes) 
  
    zipWithPad :: (A.Num a, A.Num b, Elt c) => (Exp a -> Exp b -> Exp c) -> 
        Acc (Vector a) -> Acc (Vector b) -> Acc (Vector c)
    zipWithPad f as bs =
        let I1 n = shape as
            I1 m = shape bs
            maxlen = A.max n m
            as' = as A.++ fill (index1 (maxlen-n)) 0
            bs' = bs A.++ fill (index1 (maxlen-n)) 0
        in  
        A.zipWith f as' bs'
    
     -- function to add two polynomials together
    sumVectors :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Vector Double)
    sumVectors p1 p2 = zipWithPad (+) p1 p2

    f :: Exp Double -> Exp Double
    f x = sin (cos x)

     -- Computes f (x_k) * cos (x_k)
    
    c0' :: (Exp Double -> Exp Double) -> Acc (Matrix Double) ->  Acc (Vector Double)
    c0' f nodes0 =
        let I2 _ n = shape nodes0
        in  A.map (\x -> 1.0 / (A.fromIntegral (n))*x )
        $ A.sum
        $ A.generate (shape nodes0) $ \(I2 j k) ->
            f (computeChebNode n k)
        {- --1/(1*(A.fromIntegral(n)+1))*
        (A.sum(A.map f nodes))   -}

    cj' :: (Exp Double -> Exp Double) -> Acc (Matrix Double) ->  Acc (Vector Double)
    cj' f nodesM =
        let I2 n _ = shape nodesM
        in  A.map (\x -> 2.0 / (A.fromIntegral n + 1.0)*x )
        $ A.sum
        $ A.generate (shape nodesM) $ \(I2 j k) ->
            f (computeChebNode n k)*cos (A.fromIntegral(j+1) *(2*A.fromIntegral(n)+1-2*A.fromIntegral(k))/(2*A.fromIntegral(n)+2)*pi)
         --   A.fromIntegral(j)
    
    chebCoeff' :: (Exp Double -> Exp Double) -> Exp Int -> Acc (Vector Double)
    chebCoeff' f n =
      let nodesV  = chebNodes n
          nodesM  = A.replicate (lift (Z :. n :. All)) nodesV
          nodes0  = A.replicate (lift (Z :. (1::Int) :. All)) nodesV
       in
       A.reverse ((c0' f nodes0) A.++ (cj' f nodesM ))
        

     -- Takes in an order. Returns list of chebyshev polynomials
    chebPol :: Int -> [[Double]]
    chebPol 0 = [[1.0]]
    chebPol 1 = [[0.0, 1.0], [1.0]]
    chebPol n =
        let prevResult   = chebPol (n-1)
            multTwo      = P.map (\x -> P.map (*2.0) x) prevResult
            firstTerm    = P.map (\x -> 0:x) multTwo
            subtractTerm = sumVectorsL (P.head firstTerm) (P.map P.negate (P.head (P.tail prevResult)))
        in
        subtractTerm : prevResult

    padList :: [Double] -> Int -> [Double]
    padList lst n = 
        let extraDigits = n - P.length lst
        in
            if extraDigits P.<= 0 
                then lst
                else lst P.++(P.replicate extraDigits 0)

    sumVectorsL :: [Double] -> [Double] -> [Double]
    sumVectorsL p1 p2 =
        let l1     = P.length p1
            l2     = P.length p2
            maxlen = P.max l1 l2
        in
        P.zipWith (+) (p1 P.++ P.replicate (maxlen - l1) 0)
                      (p2 P.++ P.replicate (maxlen - l2) 0)

    genChebMatrix :: Int -> Matrix Double
    genChebMatrix n = 
        let chebPolynomials = chebPol (n)  
            matrix = P.map (\x -> padList x (n+1)) chebPolynomials
            flattened = P.concat matrix
        in
        A.fromList (Z:.(n+1):.(n+1)) flattened

    multiplyCoeff :: Exp Double -> Acc (Vector Double) -> Acc (Vector Double)
    multiplyCoeff coeff vec = A.map (* coeff) vec
    
    -- Given a function f, and degree n, calculates chebyshev approximation
    -- Get list of coeffs and chebyshev polynomials. Want to zip each coeff w/ respective polynomial and multiply. 
    -- Finally, fold over all polynomials
    chebf :: (Exp Double -> Exp Double) -> Int -> Acc (Vector Double)
    chebf f n =
        let n'       = constant n
            coeffs   = chebCoeff' f n' -- size n+1 vector
            chebPols = genChebMatrix n -- size (n+1)*(n+1) matrix
            coeffsM  = A.replicate (lift (Z :. All :. (n'+1))) coeffs
        in
        A.sum $ A.transpose $ A.zipWith (*) coeffsM (use chebPols)
     --   A.replicate (lift (Z :. All :. n')) coeffs
 