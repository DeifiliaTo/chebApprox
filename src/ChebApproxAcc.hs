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
         let nodes = enumFromN (lift (Z:.n)) 0 in
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
    sumVectors p1 p2 =
        zipWithPad (+) p1 p2

    sumList :: Acc (Vector Double) -> Acc (Scalar Double)
    sumList lst = A.fold (+) 0 lst

    f :: Exp Double -> Exp Double
    f x = 2+x/(cos x)

     -- Computes f (x_k) * cos (x_k)
    computeProduct ::  (Exp Double -> Exp Double) -> Exp Int -> Exp Int -> Exp Int -> Exp Double
    computeProduct f n j k =
        f ((2*A.fromIntegral(k)+1)/(2*A.fromIntegral(n)+2)*pi) * cos (A.fromIntegral (j)*(2*A.fromIntegral(k)+1)/(2*A.fromIntegral(n)+2)*pi)
    
    c0 :: (Exp Double -> Exp Double) -> Acc (Vector Double) -> Exp Int -> Exp Double
    c0 f nodes n =
        1/(1*(A.fromIntegral(n)+1))*(the (sumList(A.map f nodes)))

    cj :: (Exp Double -> Exp Double) -> Acc (Vector Double) -> Exp Int -> Exp Int -> Exp Double
    cj f nodes n k =
        2.0/(A.fromIntegral(n)+1.0)*(the (sumList(
                    A.map (\x ->
                        f (computeChebNode n x)*cos (A.fromIntegral(k) *(2*A.fromIntegral(n)+1-2*A.fromIntegral(x))/(2*A.fromIntegral(n)+2)*pi)
                    ) (enumFromN (lift (Z:.n)) 0)
                )))
    
    chebCoeff :: (Exp Double -> Exp Double) -> Exp Int -> Acc (Vector Double)
    chebCoeff f n =
        let nodes = chebNodes n 
            enumeration = enumFromN (lift (Z:.n)) 0
        in --nodes :: Acc (Vector Double)
        --( unit((c0 f nodes n))) A.++ (A.map (\x -> cj f nodes n x)  enumeration) -- doesn't work yet! need to figure out how to scan across array.
        (A.map (\x -> cj f nodes n x) enumeration)
    

    -- Takes in an order. Returns list of chebyshev polynomials
    chebPol :: Int -> [[Double]]
    chebPol 0 = [[1.0]]
    chebPol 1 = [[0.0, 1.0], [1.0]]
    chebPol n =
        let prevResult = chebPol (n-1) in
        let multTwo = P.map (\x -> P.map (*2.0) x) prevResult in
            let firstTerm = P.map (\x -> 0:x) multTwo in
                let subtractTerm = sumVectorsL (P.head (firstTerm)) ((P.map (*(-1)) (P.head (P.tail prevResult)))) in
                    subtractTerm:(prevResult)

    padList :: [Double] -> Int -> [Double]
    padList lst n = 
        let extraDigits = n - P.length lst
        in
            if extraDigits P.<= 0 then lst
            else lst P.++(P.replicate extraDigits 0)

    sumVectorsL :: [Double] -> [Double] -> [Double]
    sumVectorsL p1 p2 =
        if (P.length p1 P.>= P.length p2)
        then P.zipWith (+) p1 (p2 P.++ P.repeat 0)
        else sumVectorsL p2 p1

    genChebMatrix :: Int -> Matrix Double
    genChebMatrix n = 
        let chebPolynomials = chebPol n 
            matrix = P.map (\x -> padList x (n+1)) chebPolynomials
            flattened = P.concat matrix
        in
            A.fromList (Z:.(n+1):.(n+1)) flattened
    
--(c0 f nodes n)