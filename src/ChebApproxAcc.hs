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

    arr2 :: Acc (Vector Double)
    arr2 = use (fromList (Z :. 3) [1, 3, 2])

    arr3 :: Acc (Vector Double)
    arr3 = use (fromList (Z :. 3) [4, 5, 1])

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
 
    -- Appends 0 to front of the list
    multiplyByX :: Acc (Vector Double)  -> Acc (Vector Double)
    multiplyByX pol = (enumFromN (lift (Z:.(1::Int))) 0) A.++ pol

 {-    genList :: 
    genList = 
        A.generate (shape nodesM) $ \(I2 j k) ->
            j+k
  -}
   
    genShiftedCoeff :: Acc (Vector Double) -> Exp Int -> Acc (Matrix Double) 
    genShiftedCoeff vec n = A.generate (index2 n (2*n)) $ \(I2 j k) -> 
        cond (j A.> k A.|| (j A.< k A.&& (k-j) A.> n))
        (constant 0)
        (vec ! (I1 (k-j)))  
    
    genCoeffMatrix :: Acc (Vector Double) -> Exp Int -> Acc (Matrix Double)
    genCoeffMatrix coeff n =
        (A.replicate (lift (Z :. All :. (2*n))) coeff)

        

    multPoly :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Vector Double)
    multPoly p1 p2 =
        let I1 n = shape p1
            I1 m = shape p1
            minDim = A.min n m
            zipped = A.zipWith (*) (genCoeffMatrix p1 minDim) (genShiftedCoeff p2 minDim)
            --p1Coeff = A.transpose (A.replicate (lift (Z :. All :. n + 1)) p1)
        in 
        A.sum $ A.transpose $ zipped
    
    arr' :: Acc (Vector Double)
    arr' = use (fromList (Z :. 5) [0..]) 

    findMax :: Acc (Vector Double) -> Exp Double
    findMax vecSlice =   (the (A.maximum vecSlice))
        --(A.fold (A.max) (vecSlice ! (I1 0)) (vecSlice))

    
    envelope :: Acc (Vector Double) -> Acc (Vector Double)
    envelope coeff =
        let mag         = A.map A.abs coeff
            I1 n        = shape coeff
            coeffMatrix = A.generate (index2 n n) $ \(I2 j k) -> 
                cond (j A.> k)
                (constant 0)
                (mag ! (I1 k))  
        in A.maximum $  coeffMatrix