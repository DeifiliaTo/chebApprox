{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE DeriveDataTypeable  #-}


module  HedgehogTest
where
  import           Hedgehog
  import qualified Hedgehog.Gen as Gen
  import qualified Hedgehog.Range as Range
  import ChebApproxAcc
  import Prelude as P
  import Data.Array.Accelerate as A
  import Data.Array.Accelerate.Debug as A
  import Data.Array.Accelerate.Interpreter as I
  import Data.Array.Accelerate.LLVM.Native as CPU
  import Data.Array.Accelerate.Array.Sugar  as Sugar( (!), Arrays, Array, Shape, Elt, DIM0, DIM1, DIM2, DIM3, Z(..), (:.)(..), fromList, size ) 
  import Data.Array.Accelerate.Test.Similar
  import Test.Tasty
  import Test.Tasty.Hedgehog
  import Data.Typeable

  




  type Run  = forall a. Arrays a => Acc a -> a
  type RunN = forall f. Afunction f => f -> AfunctionR f

  dim0 :: Gen DIM0
  dim0 = return Z

  dim1 :: Gen DIM1
  dim1 = (Z :.) <$> Gen.int (Range.linear 0 1024)

  array :: (Shape sh, Elt e) => sh -> Gen e -> Gen (Array sh e)
  array sh gen = fromList sh <$> Gen.list (Range.singleton (Sugar.size sh)) gen

  
  flt_max :: A.RealFloat a => a
  flt_max = x
    where
      n     = A.floatDigits x
      b     = A.floatRadix x
      (_,u) = A.floatRange x
      x     = A.encodeFloat (b A.^n - 1) (u - n)

  flt_min :: A.RealFloat a => a
  flt_min = x
    where
      n     = A.floatDigits x
      b     = A.floatRadix x
      (l,_) = A.floatRange x
      x     = A.encodeFloat (b A.^n - 1) (l - n - 1)

  {- except :: Gen e -> (e -> Bool) -> Gen e
  except gen f  = do
    v <- gen
    when (f v) Gen.discard
    return v

  splitEvery :: Int -> [a] -> [[a]]
  splitEvery _ [] = cycle [[]]
  splitEvery n xs =
    let (h,t) = splitAt n xs
    in  h : splitEvery n t

  splitPlaces :: [Int] -> [a] -> [[a]]
  splitPlaces []     _  = []
  splitPlaces (i:is) vs =
    let (h,t) = splitAt i vs
    in  h : splitPlaces is t -}


  prop_commutative :: Property
  prop_commutative =
    property $ do
      xs <- forAll $ (Gen.list (Range.linear 0 100) (Gen.double (Range.linearFrac 0 2)))
      ys <- forAll $ (Gen.list (Range.linear 0 100) (Gen.double (Range.linearFrac 0 2)))
      (CPU.run $ use (fromList (Z :. (100::Int)) xs)) === (CPU.run $ use (fromList (Z :. (100::Int)) xs))
      --(CPU.run $ (sumVectors ( use (fromList (Z :. 100) xs)) (use (fromList (Z :. 100 ) ys)))) === (CPU.run $ (sumVectors ( use (fromList (Z :. 100) xs)) (use (fromList (Z :. 100 ) ys))) )

  tests :: IO Bool
  tests =
    checkParallel $$(discover) 


  --test = testProperty "commutative" $ test_sqrt CPU.runN sh (e (Range.linearFrac 0 flt_max))

  test_sqrt
      :: (Shape sh, Similar e, P.Eq sh, P.Floating e, A.Floating e)
      => RunN
      -> Gen sh
      -> Gen e
      -> Property
  test_sqrt runN dim e =
    property $ do
      sh <- forAll dim
      xs <- forAll (array sh e)
      let !go = runN (A.map sqrt) in go xs ~~~ mapRef sqrt xs

  mapRef :: (Shape sh, Elt a, Elt b) => (a -> b) -> Array sh a -> Array sh b
  mapRef f xs = fromFunction (arrayShape xs) (\ix -> f (xs Sugar.! ix))

  newtype TestDouble  = TestDouble  Bool deriving (P.Eq, P.Show, Typeable)

  test_map :: RunN -> TestTree
  test_map runN =
    testGroup "map"
      [ 
        at @TestDouble $ testFloatingElt Gen.double
      ]
    where

      testFloatingElt
          :: forall a. (P.RealFloat a, A.Floating a, A.RealFrac a, Similar a)
          => (Range a -> Gen a)
          -> TestTree
      testFloatingElt e =
        testGroup (show (typeOf (undefined :: a)))
          [ 
            testDim dim1
          ]
        where
          testDim
              :: forall sh. (Shape sh, P.Eq sh)
              => Gen sh
              -> TestTree
          testDim sh =
            testGroup ("DIM" P.++ show (rank @sh))
              [ -- operators on Num
                testProperty "sqrt"       $ test_sqrt runN sh (e (Range.linearFrac 0 flt_max))
              ]

  fullrange :: P.RealFloat e => (Range e -> Gen e) -> Gen e
  fullrange gen = gen (Range.linearFracFrom 0 (-flt_max) flt_max)
      