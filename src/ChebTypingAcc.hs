{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE DeriveAnyClass  #-}
{-# LANGUAGE DeriveGeneric   #-}
{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ScopedTypeVariables #-}

module ChebTypingAcc
where
    import Data.Array.Accelerate as A
    import Data.Array.Accelerate.Debug as A
    import Data.Array.Accelerate.Interpreter as I
    import Data.Array.Accelerate.LLVM.Native as CPU
    import qualified Prelude  as P
    import ChebApproxAcc
    import qualified Data.Vector as V
    
    data Cheb = Cheb (Vector Double) deriving (Show, Generic, IsProduct Arrays, Arrays)

    x :: Acc Cheb
    x = lift (toCheb (use (fromList (Z :. 2) [0,1])))

 
    cos :: Acc Cheb -> Acc Cheb
    cos x = 
        let cosCheb = chebf (P.cos) 8 in
        lift (toCheb (composePols cosCheb (fromCheb x)))

    sin :: Acc Cheb -> Acc Cheb
    sin x = 
        let sinCheb = chebf (P.sin) 8 in
        lift (toCheb (composePols sinCheb (fromCheb x)))

    sinh :: Acc Cheb -> Acc Cheb
    sinh x = 
        let sinhCheb = chebf (P.sinh) 8 in
        lift (toCheb (composePols sinhCheb (fromCheb x)))
    
    cosh :: Acc Cheb -> Acc Cheb
    cosh x = 
        let coshCheb = chebf (P.cosh) 8 in
        lift (toCheb (composePols coshCheb (fromCheb x)))

    exp :: Acc Cheb -> Acc Cheb
    exp x = 
        let expCheb = chebf (P.exp) 8 in
        lift (toCheb (composePols expCheb (fromCheb x)))

    {- log :: Cheb -> Cheb
    log x = 
        let f = newtonApprox (P.log) 10 in
            Cheb (fnComposition f (fromCheb(x))) -}

    {- cutOffCheb :: Cheb -> Cheb
    cutOffCheb c = 
        let env = envelope (map abs c) [0..length (c)]
            plat = plateau env 2 (length env)
        in
            if plat == (length env) then Cheb (c)
            else Cheb (take plat c) -}
    
    pattern Cheb_ :: Acc (Vector Double) -> Acc Cheb
    pattern Cheb_ xs = Pattern (xs)
    {-# COMPLETE Cheb_ #-}

    fromCheb :: Acc Cheb -> Acc (Vector Double)
    fromCheb (Cheb_ x) = x

    
    toCheb :: Acc (Vector Double) -> Acc Cheb
    toCheb = Cheb_


    instance P.Num (Acc Cheb) where
        xs + ys = toCheb(sumVectors (fromCheb (xs)) (fromCheb (ys)))
        xs - ys = toCheb(sumVectors (fromCheb (xs)) (A.map (* (-1)) (fromCheb (ys))))
        xs * ys = toCheb (multPoly (fromCheb xs) (fromCheb ys))
        
        -- xs(ys)  = composePols xs ys 

    instance P.Fractional (Acc Cheb) where
        xs / ys = 
            let div       = divPol (fromCheb xs) (fromCheb ys)
                quotient  = A.afst div
                remainder = A.asnd div
                num'      = polToFn (V.fromList (CPU.run $ remainder))
                denom'    = polToFn (unlift $ (fromCheb ys))
                g (x::Exp Double)       = (num' (x::Exp Double))/(denom' (x::Exp Double))
                result    = chebfPrecise g
            in (toCheb result) + (toCheb quotient)

    instance P.Floating (Acc Cheb) where
        cos = P.id

    --formFunction :: 

    infixl 7 *^
    (*^) :: Exp Double -> Acc Cheb -> Acc Cheb
    (*^)  x ys  = toCheb (A.map (A.*x) (fromCheb ys)) 

    infixl 7 ^*
    (^*) :: Acc Cheb -> Exp Double -> Acc Cheb
    (^*) ys x = toCheb (A.map (A.*x) (fromCheb ys)) 

    infixl 7 */
    (*/) :: Exp Double -> Acc Cheb -> Acc Cheb
    (*/)  x ys  = toCheb (A.map (A./x) (fromCheb ys)) 

    infixl 7 /*
    (/*) :: Acc Cheb -> Exp Double -> Acc Cheb
    (/*) ys x = toCheb (A.map (A./x) (fromCheb ys)) 
