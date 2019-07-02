{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE DeriveAnyClass  #-}
{-# LANGUAGE DeriveGeneric   #-}
{-# LANGUAGE PatternSynonyms #-}

module ChebTypingAcc
where
    import Data.Array.Accelerate as A
    import Data.Array.Accelerate.Debug as A
    import Data.Array.Accelerate.Interpreter as I
    import Data.Array.Accelerate.LLVM.Native as CPU
    import qualified Prelude  as P
    import ChebApproxAcc
    
    data Cheb = Cheb (Acc (Vector Double)) deriving (Show, Generic)
       -- deriving (Show, Generic, IsProduct Arrays, Arrays)

    x :: Cheb
    x = Cheb (enumFromN 0 1)
    
    -- given two chebyshev functions, 
    calcCheb :: (Double -> Double) -> Cheb -> Int -> Cheb
    calcCheb f (Cheb x) n =
        let fnRep = chebf f n 
            env = envelope (map abs fnRep) [0..length (fnRep)]
        in
            if (plateau env 2 (length env)) >= (length(env) P.-1) then
                if n > 200 then Cheb (fnComposition fnRep (fromCheb x))
                else
                    calcCheb f x (n P.* 2)
            else
                cutOffCheb (fnComposition fnRep (fromCheb x))

    cos :: Cheb -> Cheb
    cos x = 
        calcCheb (P.cos) x 8

    sin :: Cheb -> Cheb
    sin x = 
        calcCheb (P.sin) x 8

    sinh :: Cheb -> Cheb
    sinh x = 
        calcCheb (P.sinh) x 8
    
    cosh :: Cheb -> Cheb
    cosh x = 
        calcCheb (P.cosh) x 8

    exp :: Cheb -> Cheb
    exp x = 
        calcCheb (P.exp) x 8

    log :: Cheb -> Cheb
    log x = 
        let f = newtonApprox (P.log) 10 in
            Cheb (fnComposition f (fromCheb(x)))

 
    --fromCheb (x::Float) = [fromx]
    --fromCheb (x::Double) = [x]
    --fromCheb (x::Int) = [fromInteger x]

    cutOffCheb :: Cheb -> Cheb
    cutOffCheb c = 
        let env = envelope (map abs c) [0..length (c)]
            plat = plateau env 2 (length env)
        in
            if plat == (length env) then Cheb (c)
            else Cheb (take plat c)
    
    pattern Cheb_ :: Acc (Vector Double) -> Acc Cheb
    pattern Cheb_ xs = Pattern (xs)
    {-# COMPLETE Cheb_ #-}

    fromCheb :: Acc Cheb -> Acc (Vector Double)
    fromCheb (Cheb_ x) = x

    toCheb :: Acc (Vector Double) -> Acc Cheb
    toCheb = Cheb_


    instance P.Num (Acc Cheb) where
        xs + ys = cutOffCheb (sumVectors xs ys)

    instance P.Fractional (Acc Cheb) -- where
        -- (/) = error "sadness"

    instance P.Floating (Acc Cheb) where
        cos = P.id

