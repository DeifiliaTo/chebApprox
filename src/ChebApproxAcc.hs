{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DeriveAnyClass  #-}
{-# LANGUAGE DeriveGeneric   #-}
{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE OverloadedStrings #-}

module  ChebApproxAcc
where
    import Prelude as P
    import Data.Array.Accelerate as A
    import Data.Array.Accelerate.Debug as A
    import Data.Array.Accelerate.Interpreter as I
    import Data.Array.Accelerate.LLVM.Native as CPU
    import qualified ChebyshevApproximations as C
    import qualified Data.Vector as V

    import Control.Exception
    import Formatting
    import Formatting.Clock
    import System.Clock

    {- arr2 :: Acc (Vector Double)
    arr2 = use (fromList (Z :. 3) [1, 3, 2]) -}

    {- value :: Exp (Int)
    value = (constant 2)

    value2 :: Exp (Int)
    value2 = (constant 1)

    value3 :: Exp (Double)
    value3 = (constant 0.02)

    arr1 :: Acc (Vector Double)
    arr1 = use (fromList (Z :. 10) [0..])


    arr3 :: Acc (Vector Double)
    arr3 = use (fromList (Z :. 3) [1, 3, 2])

    vec0 :: Acc (Vector Double)
    vec0 = use (fromList (Z :. 1) [1..])

    vec1 :: Acc (Vector Double)
    vec1 = use (fromList (Z :. 2) [0..]) -}

    -- Returns value of chebyshev zero
    computeChebNode :: Exp (Int) -> Exp (Int) -> Exp (Double)
    computeChebNode n k =
        cos ((2*A.fromIntegral(n)+1-2*A.fromIntegral(k))*pi /(2*A.fromIntegral(n)+2))

     -- Creates list of chebyshev zeros
    chebNodes :: Exp Int -> Acc (Vector Double)
    chebNodes n =
      let nodes = enumFromN (I1 (n+1)) 0
       in A.map (\x -> computeChebNode n x) nodes

    zipWithPad :: (A.Num a, A.Num b, Elt c) => (Exp a -> Exp b -> Exp c) ->
        Acc (Vector a) -> Acc (Vector b) -> Acc (Vector c)
    zipWithPad f as bs =
        let I1 n = shape as
            I1 m = shape bs
            maxlen = A.max n m
            as' = as A.++ fill (index1 (maxlen-n)) 0
            bs' = bs A.++ fill (index1 (maxlen-m)) 0
        in
        A.zipWith f as' bs'

     -- function to add two polynomials together
    sumVectors :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Vector Double)
    sumVectors p1 p2 = zipWithPad (+) p1 p2

    f :: Exp Double -> Exp Double
    f x = sin(x)

    g :: Exp Double -> Exp Double
    g x = sin(x) + sin (x*x)

     -- Computes f (x_k) * cos (x_k)

    c0' :: (Exp Double -> Exp Double) -> Acc (Matrix Double) ->  Acc (Vector Double)
    c0' f nodes0 =
        let I2 m n = shape nodes0
        in  A.map (\x -> 1.0 / (A.fromIntegral n) * x)
        $ A.sum
        $ A.generate (I2 m n) $ \(I2 j k) ->
            f (computeChebNode (n-1) k)

        {- --1/(1*(A.fromIntegral(n)+1))*
        (A.sum(A.map f nodes))   -}

    c0'L :: Exp Int
         -> (Acc (Matrix Double) -> Acc (Matrix Double))
         -> Acc (Vector Double)
    c0'L n f
      = A.map (\x -> 1.0 / (A.fromIntegral n) * x )
      $ A.sum
      $ f
      $ A.generate (I2 1 n) $ \(I2 j k) -> computeChebNode (n-1) k

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

    cj'L :: Exp Int
         -> Exp Int
         -> (Acc (Matrix Double) -> Acc (Matrix Double))
         -> Acc (Vector Double)
    cj'L n m f
      = A.map (\x -> 2.0 / (A.fromIntegral n + 1.0) * x )
      $ A.sum
      $ f
      $ A.generate (I2 n m)
      $ \(I2 j k) ->
          computeChebNode n k*cos (A.fromIntegral(j+1) *(2*A.fromIntegral(n)+1-2*A.fromIntegral(k))/(2*A.fromIntegral(n)+2)*pi)
       --   A.fromIntegral(j)


    chebCoeff' :: (Exp Double -> Exp Double) -> Exp Int -> Acc (Vector Double)
    chebCoeff' f n =
      let nodesV  = chebNodes n
          nodesM  = A.replicate (lift (Z :. n :. All)) nodesV
          nodes0  = A.replicate (lift (Z :. (1::Int) :. All)) nodesV
       in
       A.reverse ((c0' f nodes0) A.++ (cj' f nodesM ))

    chebCoeff'L
        :: (Acc (Matrix Double) -> Acc (Matrix Double))
        -> Exp Int
        -> Acc (Vector Double)
    chebCoeff'L f n =
      let nodesV = chebNodes n
          I1 m   = shape nodesV
       in
       c0'L m f A.++ cj'L n m f


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

    multiplyCoeff :: Exp Double -> Acc (Vector Double) -> Acc (Vector Double)
    multiplyCoeff coeff vec = A.map (* coeff) vec

    -- Given a function f, and degree n, calculates chebyshev approximation
    -- Get list of coeffs and chebyshev polynomials. Want to zip each coeff w/ respective polynomial and multiply.
    -- Finally, fold over all polynomials
    chebf :: (Exp Double -> Exp Double) -> Exp Int -> Acc (Vector Double)
    chebf f n =
        let coeffs   = A.reverse $ chebCoeff' f n -- size n+1 vector
            chebPols = genChebPolAcc n -- size (n+1)*(n+1) matrix
            coeffsM  = A.replicate (lift (Z :. All :. (n+1))) coeffs
        in
        A.sum $ A.transpose $ A.zipWith (*) coeffsM chebPols

    chebfL :: (Acc (Matrix Double) -> Acc (Matrix Double))
           -> Exp Int
           -> Acc (Vector Double)
    chebfL f n =
        let coeffs   = chebCoeff'L f n -- size n+1 vector
            chebPols = genChebPolAcc n -- size (n+1)*(n+1) matrix
            coeffsM  = A.replicate (lift (Z :. All :. (n+1))) coeffs
        in
        A.sum $ A.transpose $ A.zipWith (*) coeffsM chebPols

  {-
    Generates a matrix of coefficients with padded zeros. Ex: vec = [1, 2, 3]:
    [1, 2, 3, 0, 0, 0]
    [0, 1, 2, 3, 0, 0]
    [0, 0, 1, 2, 3, 0]
    [0, 0, 0, 1, 2, 3]
  -}
    genShiftedCoeff :: Acc (Vector Double) ->  Exp Int -> Acc (Matrix Double)
    genShiftedCoeff vec n  = A.generate (index2 n n) $ \(I2 j k) ->
        cond (j A.> k A.|| (j A.< k A.&& (k-j) A.>= n))
        (constant 0)
        (vec ! (I1 (k-j)))

  {-
    Generates a matrix of coefficients of same size as genShiftedCoeff. One of pols when multiplying two together.
    Ex: vec = [4, 5]:
    [4, 4, 4, 4, 4, 4]
    [5, 5, 5, 5, 5, 5]
    [0, 0, 0, 0, 0, 0]
    [0, 0, 0, 0, 0, 0]generate
  -}
    genCoeffMatrix :: Acc (Vector Double) -> Exp Int  -> Acc (Matrix Double)
    genCoeffMatrix coeff n =
        (A.replicate (lift (Z :. All :. n)) coeff)

  {-
    Multiplies two polynomials.
    Ex: p1 = [4, 5] p2 = [1, 3]:
    [4, 17, 15]
  -}
    multPoly :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Vector Double)
    multPoly p1 p2 =
        let I1 n = shape p1
            I1 m = shape p2
            dimPol = n+m-1
            zipped = A.zipWith (*) (genCoeffMatrix p1 dimPol) (genShiftedCoeff p2 dimPol)
            --p1Coeff = A.transpose (A.replicate (lift (Z :. All :. n + 1)) p1)
        in
        A.sum $ A.transpose $ zipped

    {-Example case-}
    arr' :: Acc (Vector Double)
    arr' = use (fromList (Z :. 5) [0..])

    res :: Acc (Vector Double)
    res = chebf f 25

    env :: Acc (Vector Double)
    env = envelope res

    envMat :: Acc (Matrix Double)
    envMat = genEnvMatrix env

  {-
    Given a list of coefficients, return a list of coefficients of max magnitude of descending order
    Ex: [3, -2, 0.1, 0.4] -> [3, 2, 0.4, 0.4]
  -}
    envelope :: Acc (Vector Double) -> Acc (Vector Double)
    envelope coeff =
        let mag         = A.map A.abs coeff
            I1 n        = shape coeff
            coeffMatrix = A.generate (index2 n n) $ \(I2 j k) ->
                cond (j A.> k)
                (constant 0)
                (mag ! (I1 k))
        in A.maximum $ coeffMatrix

    -- Define/set value of tolerance
    tol :: Exp Double
    tol = constant (1e-10)

  {-
    Given an envelope, generates a matrix padded with ones in bottom left corner. First two digits are cut off by definition
    Ex:
    env = [4, 2, 1, 0.1, 0.1]
    genEnvMatrix env =
        [1, 0.1, 0.1]
        [0, 0.1, 0.1]
        [0, 0.0, 0.1]
  -}
    genEnvMatrix :: Acc (Vector Double) -> Acc (Matrix Double)
    genEnvMatrix env =
        let I1 n = shape env
        in
            A.generate (index2 (n-2) (n-2)) $ \(I2 j k) ->
            cond (j A.> k)
            (constant 0)
            (env ! (I1 (k+2)))

  {-
    For all values in env, calculate 3*(1-log (x) / log (tol))
  -}
    calcR :: Acc (Vector Double) -> Acc (Vector Double)
    calcR env = A.map (\x -> 3*(1-log (x)/log (tol))) (A.drop 2 env) -- starts @ index 2 and up

    -- Example
    rl :: Acc (Vector Double)
    rl = calcR env

  {-
    Returns boolean matrix for whether or not the index of (j, k) is a valid plateau point.
    (Can probably be done more efficiently in 1 dimension)
  -}
    plateauMatrix :: Acc (Matrix Double) -> Acc (Vector Double) -> Acc (Matrix Double)
    plateauMatrix envMat rList =
        let I1 n = shape rList in
            A.generate (index2 n n) $ \(I2 j k) ->
                let kIndex = (A.round (((constant 1.25):: Exp Double )*A.fromIntegral(j)+5)) :: Exp Int
                in
                cond (k A.== kIndex A.&& (kIndex A.< n A.|| kIndex A.== n) A.&& ((envMat ! (I2 j kIndex)) / (envMat ! (I2 j 0)) A.> (rList ! (I1 j))) )
                --cond (envMat ! (I2 j k) A./= 0 A.&& j A./= 0 A.&& envMat ! (I2 j k) A.> (rList ! (I1 k)))
                (constant 1)
                (constant 0)

  {-
    Given a chebyshev function representation, return the plateau point (or size n-1 if none)
    Calculates the envelope, plateau matrix, and finds max vertically. Will either have 0 or 1 in array "summed"
    Then, zips wth an enumeration to determine the index that the plateau point is at.
  -}
    plateauPoint :: Acc (Vector Double) -> Acc (Scalar (Double))
    plateauPoint cheb =
        let envFn = envelope cheb
            I1 n = shape envFn
            envMat = genEnvMatrix envFn
            rl = calcR envFn
            pMatrix = plateauMatrix envMat rl
            summed = A.maximum $ pMatrix
            zipped = A.zipWith (*) summed (enumFromN (lift (Z:.(n+2))) 2)
            T2 filteredPos _ = (A.filter (\x -> x A.> 0) zipped)
        in
        A.minimum filteredPos

    extract :: (Vector Double, Array DIM0 Int) -> Acc (Vector Double)
    extract (x, y) = use x

    ifEnough :: Exp Double -> Acc (Vector Double) -> Acc (Scalar Bool) -- Needs to be this type signature if used in awhile
    ifEnough criteria arr =
      --let val = ( (plateauPoint arr) ! (I1 0) )in
      let pt = the (plateauPoint arr)
          I1 len = shape arr
      in
      acond
      (
        pt A.> A.fromIntegral (len)
      )
      (
        unit (constant False)
      )
      (
        acond (pt A.> criteria)
        (unit (constant True))
        (unit (constant False))
      )


    extendChebf
          :: (Exp Double -> Exp Double)
          -> Acc (Vector Double)
          -> Acc (Vector Double)
    extendChebf f chebRep =
      let I1 n = shape chebRep
       in chebf f (2*n)

    extendChebfL
        :: (Acc (Matrix Double) -> Acc (Matrix Double))
        -> Acc (Vector Double)
        -> Acc (Vector Double)
    extendChebfL f chebRep =
      let I1 n = shape chebRep
       in chebfL f (2*n)

    chebfPrecise :: (Exp Double -> Exp Double) -> Acc (Vector Double)
    chebfPrecise f =
      let orig = chebf f 5
      in
        awhile (ifEnough (constant 7)) -- TODO 30 will be changed
        (
          extendChebf f  -- By defn of awhile, need function to extend (Type Acc (Vector Double) -> Acc (Vector Double))
        )
        (
          orig
        )

    chebfPreciseL
        :: (Acc (Matrix Double) -> Acc (Matrix Double))
        -> Acc (Vector Double)
    chebfPreciseL f =
      let orig = chebfL f 5
      in
        awhile (ifEnough (constant 7)) -- TODO 30 will be changed
        (
          extendChebfL f  -- By defn of awhile, need function to extend (Type Acc (Vector Double) -> Acc (Vector Double))
        )
        (
          orig
        )

    -- given a matrix, compute an additional row
    -- Take in the matrix computed so far, total dim of matrix, and row we are currently on
    chebPolAcc :: Exp Int -> Acc (Matrix Double) -> Acc (Matrix Double)
    chebPolAcc n mat =
      let I2 j _  = shape mat
          nextRow = generate (index2 1 n) $ \(I2 one k) ->
            cond (k A.> j)
            (
              0
            )
            (
              cond (j A.== 0)
              (
                -1*(mat ! index2 (j-2) k)
              )
              (
                2*(mat ! index2 (j-1) (k-1)) - (mat ! index2 (j-2) k)
              )
            )
          nextRow' = A.transpose nextRow
          mat' = A.transpose mat
      in
      A.transpose $ (mat' A.++ nextRow')

    chebPolAccBase :: Exp Int -> Acc (Matrix Double)
    chebPolAccBase n =
      generate (index2 2 n) $ \(I2 j k) ->
        cond (j A.== k)
        (1)
        (0)


    ifIter :: Exp Int -> Acc (Matrix Double) -> Acc (Scalar Bool)
    ifIter n mat =
      let I2 j _  = shape mat in
      acond (j A.< n)
      (unit (constant True))
      (unit (constant False))

    genChebPolAcc :: Exp Int -> Acc (Matrix Double)
    genChebPolAcc n =
      let n'   = n + 1
          base = chebPolAccBase n'
      in
        awhile (ifIter n')
        (chebPolAcc n')
        (base)

    samplePol :: Acc (Vector Double)
    samplePol = use (fromList (Z :. 3) [1, 0, 2]) -- 1 + x ^2
    -- if coeffs are: [1, 0, 2], I need to get vals for [-1, 0, 1, 2], then:
    -- [1, 0, 2] [-1 -1 -1]
    -- [1, 0, 2] [ 0  0  0]
    -- [1, 0, 2] [ 1  1  1]
    -- [1, 0, 2] [ 2  2  2]
    evalOverRange :: Acc (Vector Double) -> Exp Int -> Acc (Vector Double)
    evalOverRange pol numPoints =
      let I1 n     = shape pol
          lst      = enumFromStepN (lift (Z :. (numPoints))) (constant (-1)) (2.0/(A.fromIntegral(numPoints)-1))
          repCoeff =  A.transpose $ A.replicate (lift (Z :. n :. All)) lst
          eval     = A.generate (index2 (numPoints) n) $ \(I2 j k) ->
            pol ! (I1 k) *
            (repCoeff ! (I2 j k) A.^ k)
      in
        A.sum $ eval

    evalFunction :: (Exp Double -> Exp Double) -> Exp Int -> Acc (Vector Double)
    evalFunction f numPoints =
      let lst = enumFromStepN (lift (Z :. (numPoints))) (constant (-1)) (2.0/(A.fromIntegral(numPoints)-1))
      in
        A.map f lst

    approxError :: (Exp Double -> Exp Double) -> Acc (Vector Double) -> Acc (Vector Double)
    approxError f pol =
      let numPoints = constant (100::Int)
          approx    = evalOverRange pol numPoints
          real      = evalFunction f numPoints
          width     = constant 2.0/((A.fromIntegral(numPoints)::Exp Double)-1.0)
      in
      A.map (A.* width) (A.map (A.^(2::Exp Int))(A.map (A.abs) (A.zipWith (-) real approx )))

    composeMat1 :: Acc (Vector Double) -> Exp Int -> Acc (Matrix Double)
    composeMat1 pol1 m =
      let I1 n = shape pol1
      in
        A.transpose $ A.replicate (lift (Z :. (n*m-2) :. All)) pol1


    ifIterPow :: Exp Int -> Acc (Vector Double) -> Acc (Scalar Bool)
    ifIterPow totalSize vec =
      let I1 j = shape vec in
      acond (j A.< totalSize)
      (unit (constant True))
      (unit (constant False))

    -- Multiply something of size n, m times.
    pow :: Exp Int -> Exp Int -> Acc (Vector Double) -> Acc (Vector Double)
    pow nSize mTimes vec =
      let totalSize = (nSize-1) A.* mTimes-2 in
      awhile (ifIterPow totalSize)
      (multPoly vec)
      (vec)

   -- given a matrix, compute an additional row
    -- Take in the matrix computed so far, total dim of matrix, and row we are currently on
    chebCompAcc :: Exp Int -> Acc (Vector Double) -> Acc (Matrix Double) -> Acc (Matrix Double)
    chebCompAcc n vec mat =
      let I2 m _  = shape mat
          nextRow = pow n (m) vec
          nextRow' = A.transpose $ A.generate (index2 1 ((n-1)*m+1)) $ \(I2 j k) ->
            cond (k A.>= (n-1)*(m)-1) -- TODO check dimensions
            (0)
            (nextRow ! (I1 k))
          mat' = A.transpose mat
      in
      A.transpose $ (mat' A.++ nextRow')

    {-
        Generates the first two rows (padded with zeros to the appropriate length) of Matrix req for
        function composition
    -}
    chebCompAccBase :: Exp Int -> Exp Int -> Acc (Vector Double) -> Acc (Matrix Double)
    chebCompAccBase n m vec =
      A.generate (index2 2 ((n-1)*m+1)) $ \(I2 j k) ->
        let I1 m = shape vec in
          cond (j A.== 0)
          (
            cond (k A./= 0)
            (
              0
            )
            (
              1
            )
          )
          (
            cond (k A.< m)
            (
              vec ! (I1 k)
            )
            (0)
          )

    {-
          If statement to see if it needs to keep generating rows
    -}
    ifIterComp :: Exp Int -> Acc (Matrix Double) -> Acc (Scalar Bool)
    ifIterComp n mat =
      let I2 j _  = shape mat in
      acond (j A.< n)
      (unit (constant True))
      (unit (constant False))

    {-
        Generate matrix2 required for composition
    -}
    genChebCompAcc :: Exp Int -> Exp Int -> Acc (Vector Double) -> Acc (Matrix Double)
    genChebCompAcc n m vec =
      let n'   = n + 1
          m'   = m + 1
          base = chebCompAccBase n' m' vec
      in
        awhile (ifIterComp n')
        (chebCompAcc n' vec)
        (base)

    {-
        Function for polynomial composition
    -}
    composePols :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Vector Double)
    composePols p1 p2 =
      let I1 n = shape p1
          I1 m = shape p2
          m1   = composeMat1 p1 m
          m2   = genChebCompAcc n m p2
      in
      A.sum $ A.transpose $ A.zipWith (*) m1 m2

    dPol' :: Acc (Vector Double)
    dPol' = use (fromList (Z :. 5) [12, 32, 27, 11, 2])

    nPol' :: Acc (Vector Double)
    nPol' = use (fromList (Z :. 3) [2, 5, 2])

    eps :: Exp Double
    eps = constant 1e-10

    genDivBase :: Exp Int -> Exp Int -> Acc (Vector Double) -> Acc (Vector Double) -> Acc (Matrix Double)
    genDivBase d n dPol nPol =
      A.generate (index2 2 (d+1)) $ \(I2 t j) ->
        A.cond (t A.== 0)
        (
          A.cond (A.fromIntegral(j) A.<= (A.fromIntegral (n)+eps))
          (
            nPol ! (I1 (n-j))
          )
          (
            0
          )
        )
        (
          dPol ! (I1 (d-j))
        )

    genDivMatrix :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Matrix Double)
    genDivMatrix dPol nPol =
      let an   = nPol ! (I1 0)
          I1 n = shape nPol   -- 3
          I1 d = shape dPol   -- 5
          base = genDivBase (d-1) (n-1) dPol' nPol'
      in
        awhile (ifIterComp (d-n+3))
        (genDivRow (d-1) (n-1) an nPol)
        (base)

    genRemainder :: Acc (Matrix Double) -> Exp Int -> Acc (Vector Double)
    genRemainder mat len =
      let I2 rows _ = shape mat
      in
        A.generate (index1 len) $ \ (I1 j) ->
          mat ! (I2 (rows-1) j)

    genQuotient :: Acc (Matrix Double) -> Exp Int -> Acc (Vector Double)
    genQuotient mat len =
      let I2 rows _ = shape mat
          an        = mat ! (I2 0 0)
      in
        A.generate (index1 len) $ \ (I1 j) ->
          mat ! (I2 (j+1) 0) / an

    divPol :: Acc (Vector Double) -> Acc (Vector Double) -> Acc ((Vector Double), (Vector Double))
    divPol dPol nPol =
      let I1 n      = shape nPol
          I1 d      = shape dPol
          mat       = genDivMatrix dPol nPol
          remainder = genRemainder mat (n-1) -- n-1 = degree of pol
          quotient  = genQuotient mat (d-n+1)
      in
      lift (quotient, remainder)


    -- given a matrix, compute an additional row
    -- Take in the matrix computed so far, total dim of matrix, and row we are currently on
    genDivRow :: Exp Int -> Exp Int -> Exp Double -> Acc (Vector Double) -> Acc (Matrix Double) -> Acc (Matrix Double)
    genDivRow d n an nPol mat =
      let
        I2 currRow _ = shape mat -- 2 x something
        nextRow =  A.transpose $ A.generate (index2 1 (d+1)) $ \(I2 t j) ->
          A.cond ((n-1-j) A.>= 0)
          (
            (an * mat ! (I2 (currRow -1) (j+1)) - (nPol ! (I1 (n-1-j)))*(mat ! (I2 (currRow-1) 0)))/an
          )
          (
            --A.cond ((j - currRow +1 )A.< d)
            A.cond ((d-j) A.< currRow-1)

            (0)
            (
              (an * mat ! (I2 (currRow-1) (j+1)))/an
            )
          )
      in
      A.transpose $ ((A.transpose $ mat) A.++ nextRow)

    clockSomething :: a -> IO ()
    clockSomething something =
      do
        start <- getTime Monotonic
        evaluate (something)
        end <- getTime Monotonic
        fprint (timeSpecs % "\n") start end

    {-
        Takes a polynomial and converts it to a function.
    -}
  {-   polToFn :: (V.Vector Double) -> (Exp Double -> Exp Double)
    polToFn pol =
      let f' (x::Exp Double) = evalPol' pol x
      in
        f'

    computeVal' :: Double -> Double -> Int -> Double
    computeVal' x coeff rep =
      coeff * (x P.^ rep)

    -- given pol and point, compute
    evalPol' :: (V.Vector Double) -> Double -> Exp Double
    evalPol' pol x =
      let n      = V.length pol
          zipped = V.zipWith (computeVal' x) pol (V.enumFromN (0) (n+1)) -- 0 is starting point, n+1 is number of digits. Check to see if should be only up to n.
      in
         constant (P.foldl (+) (0) zipped) -- (Exp Double)

    fn' x = 1 + 2*x + 4*x*x
    -}
    polToFn :: Acc (Vector Double) -> (Exp Double -> Acc (Scalar Double))--(Exp Double -> Exp Double)
    polToFn = evalPol
      -- let f' (x::Exp Double) = evalPol pol x
      -- in
      --   f'

    --computeVal :: Exp Double -> Exp Double -> Exp Int -> Exp Double
    computeVal :: Exp Double -> Exp Double -> Exp Int -> Exp Double
    computeVal x coeff rep =
      coeff * (x A.^ rep)

    -- given pol and point, compute
    evalPol :: Acc (Vector Double) -> Exp Double -> Acc (Scalar Double)
    evalPol pol x =
      let I1 n   = shape pol
          zipped = A.zipWith (computeVal x) pol (enumFromN (I1 (n+1)) 0)
      in
      {- the -} A.fold (+) (0) zipped -- (Exp Double)

    evalPolL
        :: Acc (Vector Double)  -- current polynomial coefficient
        -> Acc (Matrix Double)  -- new points to evaluate at
        -> Acc (Matrix Double)  -- new computed values
    evalPolL pol xs =
      let I1 c   = shape pol
          I2 n m = shape xs
          pol'   = A.replicate (lift (Z :. n :. m :. All)) pol
          xs'    = A.replicate (lift (Z :. All :. All :. c)) xs
          zipped = A.izipWith (\(I3 _ _ i) x p -> computeVal x p i) xs' pol'
       in
       A.sum zipped

