module ChebyshevApproximations 
where
    import Prelude
    
    -- Define data type to be able to build complex expressions
    data Expr =  Value Float | Cos Expr | Sin Expr | Tan Expr | Exp Expr | Log Expr 

    -- TODO: Change type to floating (more general)
    {-
    
    -}

    -- sin (cos (x)) = eval (Compute Sin (Compute Cos (Value x)))
    eval :: Expr -> Float
    eval (Value x) = x
    eval (Cos x) = cos (eval (x))
    eval (Sin x) = sin (eval (x))
    eval (Tan x) = sin (eval x) / cos (eval x)
    eval (Exp x) = exp (eval (x))
    eval (Log x) = log (eval (x))
    
    -- 0 : Newtonian approach
    -- 1 : lagrangian approach
    polConst :: Expr -> Float
    polConst (Value x) = 0
    polConst (Sin x) = 0
    polConst (Log x) = 1
    
    -- Returns value of chebyshev zero
    computeChebNode :: Float -> Float -> Float
    computeChebNode n k =
         cos ((2*n+1-2*k)*pi/(2*n+2))

    -- Creates list of chebyshev zeros
    chebNodes :: Float  -> [Float]
    chebNodes n = 
        let nodes = [0..n] in
            map (\x -> computeChebNode n x) (nodes)
    
    finalCheb :: [Float] -> [Float] -> [Float]
    finalCheb (x:xs) (y:ys) =
        if xs == [] then (x-y):ys
        else if ys  == [] then (x-y):xs
        else (x-y):(finalCheb xs ys)

    -- Takes in an order. Returns list of chebyshev polynomials
    chebPol :: Float  -> [[Float]]
    chebPol 0 = [[1]]
    chebPol 1 = [[0, 1], [1]]
    chebPol n = 
        let prevResult = chebPol (n-1) in
        let multTwo = map (\x -> map (*2) x) prevResult in
            let firstTerm = map (\x -> 0:x) multTwo in
                let subtractTerm = finalCheb (head (firstTerm)) (head (tail prevResult)) in
                  subtractTerm:(prevResult)

    
    -- calcCoeffs :: ([Float] -> [Float]) -> Float -> Float
    -- will be much easier w/ vector operation
    -- Adds up all the elements in a vector
    -- TODO: Replace with vector ops/accelerate library
    sumList :: [Float] -> Float
    sumList lst = foldl (+) 0 lst

    -- Computes f (x_k) * cos (x_k)
    computeProduct :: (Float -> Float) -> Float -> Float -> Float -> Float
    computeProduct f n j k = 
        f ((2*k+1)/(2*n+2)*pi) * cos (j*(2*k+1)/(2*n+2)*pi)

    -- Calculates c_j (coefficient)
    chebCoeff :: (Float -> Float) -> Float  -> Float -> [Float]
    chebCoeff f n iter = 
        let nodes = chebNodes n in
            if iter == 0 then [1/(1*(n+1))*sumList (map f nodes)]
            else 
                2/(n+1)*sumList(
                    map (\x -> 
                        --computeProduct f n iter x
                        f (computeChebNode n x)*cos (iter *(2*n+1-2*x)/(2*n+2)*pi)
                    ) [0..n]
                ):
                (chebCoeff f n (iter-1))
    
    -- function to add two polynomials together
    sumVectors :: [Float] -> [Float] -> [Float]
    sumVectors p1 p2 = 
        if (length p1 >= length p2)
        then zipWith (+) p1 (p2 ++ repeat 0)
        else sumVectors p2 p1
    

    -- takes in list of coefficients, and cheb polynomials, outputs polynomial approximation
    -- Used for c_j * T_j
    sumPoly :: [Float] -> [[Float]] -> [Float]
    sumPoly (c:cs) (p:pols) = 
        if cs == [] then
            map (*c) p
        else
            sumVectors (map (*c) p) (sumPoly cs pols)


    chebf :: (Float -> Float) -> Float -> [Float]
    chebf f n = 
        let coeffs = chebCoeff f n n in
            let chebPols = chebPol n in
                sumPoly coeffs chebPols

    polCalc ::[Float] -> Float -> Int -> Float
    polCalc (c:coeffs) x order =
        if coeffs == [] then
            c * x^order
        else
            c * x^order + polCalc coeffs x (order + 1)
    
   
    
    -- Example functions
    f :: Float -> Float
    f x = cos x

    h :: Float -> Float
    h x = sin x + exp (cos x) - x*x*x
    
    
    -- Newtonian approach
    

    findDiff :: Float -> Float -> Float -> Float
    findDiff a b order = (b-a)/order

    -- Given [1, 2, 4, 8], order 1 -> [1, 2, 4] 
    divDiff :: [Float] -> Float -> [Float]
    divDiff (x:xs) order =
        if xs == [] then []
        else 
            let b = head xs in
                (findDiff x b order):(divDiff xs order)

    
    -- Given [1, 2, 4, 8] --> [[1, 2, 4, 8], [1, 2, 4], [0.5, 1], [1/6]]
    divDiffList :: [Float] -> Float -> [[Float]]
    divDiffList (x:xs) order =
       if length (x:xs) == 2 then [divDiff (x:xs) (order)]
       else 
            let lst = (divDiff (x:xs) order) in
                lst:(divDiffList lst (order+1)) 
        
    newtonCoeffs :: [Float] -> [Float]
    newtonCoeffs lst =
        (head lst):(map head (divDiffList lst 1))

    multiplyByX :: [Float] -> [Float]
    multiplyByX p = 0:p
    
    multPoly :: [Float] -> [Float] -> [Float]
    multPoly [] p2 = []
    multPoly (p:p1) p2 = 
        let pTimesP2 = map (*p) p2 in
            let xTimesP1Timesp2 = multiplyByX (multPoly p1 p2) in
                sumVectors pTimesP2 xTimesP1Timesp2 
    
    -- Generate list of polynomials to be multiplied (Newtonian form)
    newtonPolMult :: [Float] -> [[Float]]
    newtonPolMult nodes = 
        map (\x -> [-1*x, 1]) nodes
    
    polMultFoldable :: [[Float]] -> [[Float]]
    polMultFoldable lst = 
        [1]:scanl1 (\x y -> multPoly x y) lst
    
    fn :: Float -> Float
    fn x = log x

    newtonApprox :: (Float -> Float) -> Float -> [Float]
    newtonApprox f n =
        let newtonNodes = [0.1,0.2..1]
            coeffs = newtonCoeffs (map f newtonNodes)
            pairedList = zip coeffs (polMultFoldable (newtonPolMult newtonNodes))
            mapped = map (\(x, y) -> map (*x) y) pairedList 
            in foldl (\x y -> sumVectors x y) [] mapped

    printErrLag :: (Float -> Float) -> Float -> [Float]
    printErrLag f order =
        let soln = chebf f order in
            map (\x -> (polCalc soln x 0) - f x ) [-1, -0.9..1]

    printErrNew :: (Float -> Float) -> Float -> [Float]
    printErrNew f order =
        let soln = newtonApprox f order in
            map (\x -> (polCalc soln x 0) - f x ) [0.1, 0.15..1]

