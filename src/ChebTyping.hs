module ChebTyping 
where
    import Prelude 
    import ChebyshevApproximations
    
    
    data Cheb = Cheb [Double] deriving (Show)
    
    x :: Cheb
    x = Cheb [0.0, 1.0]
    
    -- given two chebyshev functions, 
    calcCheb :: (Double -> Double) -> Cheb -> Int -> Cheb
    calcCheb f x n =
        let fnRep = chebf f n 
            env = envelope (map abs fnRep) [0..length (fnRep)]
        in
            if (plateau env 2 (length env)) >= (length(env)-1) then
                if n > 100 then Cheb (fnComposition fnRep (fromCheb x))
                else
                    calcCheb f x (n*2)
            else
                cutOffCheb (fnComposition fnRep (fromCheb x))

    cos :: Cheb -> Cheb
    cos x = 
        calcCheb (Prelude.cos) x 8

    sin :: Cheb -> Cheb
    sin x = 
        calcCheb (Prelude.sin) x 8

    sinh :: Cheb -> Cheb
    sinh x = 
        calcCheb (Prelude.sinh) x 8
    
    cosh :: Cheb -> Cheb
    cosh x = 
        calcCheb (Prelude.cosh) x 8

    exp :: Cheb -> Cheb
    exp x = 
        calcCheb (Prelude.exp) x 8

    log :: Cheb -> Cheb
    log x = 
        let f = newtonApprox (Prelude.log) 10 in
            Cheb (fnComposition f (fromCheb(x)))

    fromCheb :: Cheb -> [Double]
    fromCheb (Cheb x) = x

    cutOffCheb :: [Double] -> Cheb
    cutOffCheb c = 
        let env = envelope (map abs c) [0..length (c)]
            plat = plateau env 2 (length env)
        in
            if plat == (length env) then Cheb (c)
            else Cheb (take plat c)

    instance (Num Cheb) where
        (Cheb x) + (Cheb y) = 
            let sum = (sumVectors x y) in
                cutOffCheb sum
        (Cheb x) - (Cheb y) = 
            let diff = (sumVectors x (map (*(-1)) y)) in
                cutOffCheb diff
        (Cheb x) * (Cheb y) = 
            let mult = multPoly x y in
                cutOffCheb mult
