module ChebTyping 
where
    import Prelude 
    import ChebyshevApproximations
    
    
    data Cheb = Cheb [Double] deriving (Show)
    
    x :: Cheb
    x = Cheb [0.0, 1.0]
    
    cos :: Cheb -> Cheb
    cos x = 
        let f = chebf (Prelude.cos) 10 in
            Cheb (fnComposition f (fromCheb(x)))

    sin :: Cheb -> Cheb
    sin x = 
        let f = chebf (Prelude.sin) 10 in
            Cheb (fnComposition f (fromCheb(x)))

    sinh :: Cheb -> Cheb
    sinh x = 
        let f = chebf (Prelude.sinh) 10 in
            Cheb (fnComposition f (fromCheb(x)))
    
    cosh :: Cheb -> Cheb
    cosh x = 
        let f = chebf (Prelude.cosh) 10 in
            Cheb (fnComposition f (fromCheb(x)))

    exp :: Cheb -> Cheb
    exp x = 
        let f = chebf (Prelude.exp) 10 in
            Cheb (fnComposition f (fromCheb(x)))

    log :: Cheb -> Cheb
    log x = 
        let f = newtonApprox (Prelude.log) 10 in
            Cheb (fnComposition f (fromCheb(x)))

    fromCheb :: Cheb -> [Double]
    fromCheb (Cheb x) = x

    instance (Num Cheb) where
        (Cheb x) + (Cheb y) = Cheb (sumVectors x y)
        (Cheb x) * (Cheb y) = Cheb (multPoly x y)
    
