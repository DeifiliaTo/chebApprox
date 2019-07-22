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
    
    data Cheb = Cheb (Vector Double) deriving (Show, Generic, IsProduct Arrays, Arrays)

    x :: Acc Cheb
    x = lift (toCheb (use (fromList (Z :. 2) [0,1])))

 
    cos :: Acc Cheb -> Acc Cheb
    cos x = 
        let cosCheb = chebf (P.cos) 8 in
        lift (toCheb (composePols cosCheb (fromCheb x)))

    sin :: Acc Cheb -> Acc Cheb
    sin x = 
        lift (toCheb (chebf (P.sin) 8))

    sinh :: Acc Cheb -> Acc Cheb
    sinh x = 
        lift (toCheb (chebf (P.sinh) 8))
    
    cosh :: Acc Cheb -> Acc Cheb
    cosh x = 
        lift (toCheb (chebf (P.cosh) 8))

    exp :: Acc Cheb -> Acc Cheb
    exp x = 
        lift (toCheb (chebf (P.exp) 8))

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
        -- xs(ys)  = composePols xs ys 

    instance P.Fractional (Acc Cheb) -- where
        -- (/) = error "sadness"

    instance P.Floating (Acc Cheb) where
        cos = P.id

