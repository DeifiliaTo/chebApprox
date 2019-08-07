{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE DeriveDataTypeable  #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}


module HedgehogTest where

import ChebApproxAcc
import ChebTypingAcc

import           Hedgehog
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range
import Prelude as P
import Data.Array.Accelerate as A
import Data.Array.Accelerate.Debug as A
import Data.Array.Accelerate.Interpreter as I
import Data.Array.Accelerate.LLVM.Native as CPU
import Data.Array.Accelerate.Array.Sugar as Sugar( (!), Arrays, Array, Shape, Elt, DIM0, DIM1, DIM2, DIM3, Z(..), (:.)(..), fromList, size )
import Data.Array.Accelerate.Test.Similar
import Test.Tasty
import Test.Tasty.Hedgehog
import Data.Typeable


{--
data Expr =
    Var String
  | Lam String Expr
  | App Expr Expr

-- Assuming we have a name generator
genName :: MonadGen m => m String

-- We can write a generator for expressions
genExpr :: MonadGen m => m Expr
genExpr =
  Gen.recursive Gen.choice [
      -- non-recursive generators
      Var <$> genName
    ] [
      -- recursive generators
      Gen.subtermM genExpr (x -> Lam <$> genName <*> pure x)
    , Gen.subterm2 genExpr genExpr App
    ]
--}


-- Operator terms we want to generate random combinations of
data Term
  = Val
  | Sin Term
  | Cos Term
  | Exponential Term
  | Divide Term Term
  deriving Show

-- we can write a generator for expressions
genTerm :: MonadGen m => m Term
genTerm =
  Gen.recursive Gen.choice
    -- non-recursive generators
    [ pure Val ]
    -- recursive generators
    [ Gen.subtermM genTerm (pure . Sin)
    , Gen.subtermM genTerm (pure . Cos)
    , Gen.subtermM genTerm (pure . Exponential)
    ]

evalTerm :: Term -> Double -> Double
evalTerm Val x             = x
evalTerm (Sin t) x         = P.sin (evalTerm t x)
evalTerm (Cos t) x         = P.cos (evalTerm t x)
evalTerm (Exponential t) x = P.exp (evalTerm t x)

termToAcc :: Term -> Exp Double -> Exp Double
termToAcc Val x             = x
termToAcc (Sin t) x         = A.sin (termToAcc t x)
termToAcc (Cos t) x         = A.cos (termToAcc t x)
termToAcc (Exponential t) x = A.exp (termToAcc t x)


type Run  = forall a. Arrays a => Acc a -> a
type RunN = forall f. Afunction f => f -> AfunctionR f

dim0 :: Gen DIM0
dim0 = return Z

dim1 :: Gen DIM1
dim1 = (Z :.) <$> Gen.int (Range.linear 0 1024)

f64 :: Gen Double
f64 = Gen.double (Range.linearFracFrom 0 (-log_flt_max) log_flt_max)

array :: (Shape sh, Elt e) => sh -> Gen e -> Gen (Array sh e)
array sh gen = fromList sh <$> Gen.list (Range.singleton (Sugar.size sh)) gen

log_flt_max :: P.RealFloat a => a
log_flt_max = log flt_max

flt_max :: P.RealFloat a => a
flt_max = x
  where
    n     = P.floatDigits x
    b     = P.floatRadix x
    (_,u) = P.floatRange x
    x     = P.encodeFloat (b P.^n - 1) (u - n)

flt_min :: P.RealFloat a => a
flt_min = x
  where
    n     = P.floatDigits x
    b     = P.floatRadix x
    (l,_) = P.floatRange x
    x     = P.encodeFloat (b P.^n - 1) (l - n - 1)

-- except :: Gen e -> (e -> Bool) -> Gen e
-- except gen f  = do
--   v <- gen
--   when (f v) Gen.discard
--   return v

prop_commutative :: Property
prop_commutative =
  property $ do
    sh <- forAll dim1
    xs <- forAll (array sh f64)
    ys <- forAll (array sh f64)
    let go = CPU.runN (\x y -> fromCheb (toCheb x + toCheb y))
    go xs ys === go ys xs

prop_approx_error :: Property
prop_approx_error =
  property $ do
    term <- forAll genTerm
    let f   = termToAcc term
        pol = chebf f 20 -- TODO
        -- pol = chebfPrecise f -- TODO
        err = approxError f pol
        ok  = A.all (A.< 1.0E-15) err
    --
    True === indexArray (CPU.run ok) Z

fullrange :: P.RealFloat e => (Range e -> Gen e) -> Gen e
fullrange gen = gen (Range.linearFracFrom 0 (-flt_max) flt_max)

tests :: IO Bool
tests =
  checkParallel $$(discover)

