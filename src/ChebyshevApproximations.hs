module ChebyshevApproximations
where
    import Prelude

    -- Define data type to be able to build complex expressions
    data Expr =  Value Double | Cos Expr | Sin Expr | Tan Expr | Exp Expr | Log Expr

    -- TODO: Change type to Doubleing (more general)
    {-

    -}

    -- sin (cos (x)) = eval (Compute Sin (Compute Cos (Value x)))
    eval :: Expr -> Double
    eval (Value x) = x
    eval (Cos x) = cos (eval (x))
    eval (Sin x) = sin (eval (x))
    eval (Tan x) = sin (eval x) / cos (eval x)
    eval (Exp x) = exp (eval (x))
    eval (Log x) = log (eval (x))

    -- 0 : Newtonian approach
    -- 1 : lagrangian approach
    polConst :: Expr -> Double
    polConst (Value x) = 0
    polConst (Sin x) = 0
    polConst (Log x) = 1

    -- Returns value of chebyshev zero
    computeChebNode :: Double -> Double -> Double
    computeChebNode n k =
         cos ((2*n+1-2*k)*pi/(2*n+2))

    -- Creates list of chebyshev zeros
    chebNodes :: Double  -> [Double]
    chebNodes n =
        let nodes = [0..n] in
            map (\x -> computeChebNode n x) (nodes)

    finalCheb :: [Double] -> [Double] -> [Double]
    finalCheb (x:xs) (y:ys) =
        if xs == [] then (x-y):ys
        else if ys  == [] then (x-y):xs
        else (x-y):(finalCheb xs ys)

    -- Takes in an order. Returns list of chebyshev polynomials
    chebPol :: Double  -> [[Double]]
    chebPol 0 = [[1]]
    chebPol 1 = [[0, 1], [1]]
    chebPol n =
        let prevResult = chebPol (n-1) in
        let multTwo = map (\x -> map (*2) x) prevResult in
            let firstTerm = map (\x -> 0:x) multTwo in
                let subtractTerm = finalCheb (head (firstTerm)) (head (tail prevResult)) in
                  subtractTerm:(prevResult)


    -- calcCoeffs :: ([Double] -> [Double]) -> Double -> Double
    -- will be much easier w/ vector operation
    -- Adds up all the elements in a vector
    -- TODO: Replace with vector ops/accelerate library
    sumList :: [Double] -> Double
    sumList lst = foldl (+) 0 lst

    -- Computes f (x_k) * cos (x_k)
    computeProduct :: (Double -> Double) -> Double -> Double -> Double -> Double
    computeProduct f n j k =
        f ((2*k+1)/(2*n+2)*pi) * cos (j*(2*k+1)/(2*n+2)*pi)

    -- Calculates c_j (coefficient)
    chebCoeff :: (Double -> Double) -> Double  -> Double -> [Double]
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
    sumVectors :: [Double] -> [Double] -> [Double]
    sumVectors p1 p2 =
        if (length p1 >= length p2)
        then zipWith (+) p1 (p2 ++ repeat 0)
        else sumVectors p2 p1


    -- takes in list of coefficients, and cheb polynomials, outputs polynomial approximation
    -- Used for c_j * T_j


    chebf :: (Double -> Double) -> Double -> [Double]
    chebf f n =
        let coeffs = chebCoeff f n n
            chebPols = chebPol n
            zipped = zip coeffs chebPols
            mapped = map (\(x, y) -> map (*x) y) zipped
            in foldl (\x y -> sumVectors x y) [] mapped

    polCalc ::[Double] -> Double -> Int -> Double
    polCalc (c:coeffs) x order =
        if coeffs == [] then
            c * x^order
        else
            c * x^order + polCalc coeffs x (order + 1)



    -- Example functions


    h :: Double -> Double
    h x = sin x + exp (cos x) - x*x*x


    -- Newtonian approach
    findDiff :: Double -> Double -> Double -> Double
    findDiff a b order = (b-a)/order

    -- Given [1, 2, 4, 8], order 1 -> [1, 2, 4]
    divDiff :: [Double] -> Double -> [Double]
    divDiff (x:xs) order =
        if xs == [] then []
        else
            let b = head xs in
                (findDiff x b order):(divDiff xs order)

    -- Given [1, 2, 4, 8] --> [[1, 2, 4, 8], [1, 2, 4], [0.5, 1], [1/6]]
    divDiffList :: [Double] -> Double -> [[Double]]
    divDiffList (x:xs) order =
       if length (x:xs) == 2 then [divDiff (x:xs) (order)]
       else
            let lst = (divDiff (x:xs) order) in
                lst:(divDiffList lst (order+1))

    newtonCoeffs :: [Double] -> [Double]
    newtonCoeffs lst =
        (head lst):(map head (divDiffList lst 1))

    multiplyByX :: [Double] -> [Double]
    multiplyByX p = 0:p

    multPoly :: [Double] -> [Double] -> [Double]
    multPoly [] p2 = []
    multPoly (p:p1) p2 =
        let pTimesP2 = map (*p) p2 in
            let xTimesP1Timesp2 = multiplyByX (multPoly p1 p2) in
                sumVectors pTimesP2 xTimesP1Timesp2

    -- Generate list of polynomials to be multiplied (Newtonian form)
    newtonPolMult :: [Double] -> [[Double]]
    newtonPolMult nodes =
        map (\x -> [-1*x, 1]) nodes

    polMultFoldable :: [[Double]] -> [[Double]]
    polMultFoldable lst =
        [1]:scanl1 (\x y -> multPoly x y) lst

    fn :: Double -> Double
    fn x = 2**x

    newtonApprox :: (Double -> Double) -> Double -> [Double]
    newtonApprox f n =
        --let newtonNodes = chebNodes n
        let newtonNodes = [0.5, 0.6..2]
            coeffs = newtonCoeffs (map f newtonNodes)
            pairedList = zip coeffs (polMultFoldable (newtonPolMult newtonNodes))
            mapped = map (\(x, y) -> map (*x) y) pairedList
            in foldl (\x y -> sumVectors x y) [] mapped

    printErrLag :: (Double -> Double) -> Double -> [Double]
    printErrLag f order =
        let soln = chebf f order in
            map (\x -> (polCalc soln x 0) - f x ) [-1, -0.9..1]

    printErrNew :: (Double -> Double) -> Double -> [Double]
    printErrNew f order =
        let soln = newtonApprox f order in
            map (\x -> (polCalc soln x 0) - f x ) [0.5, 1..2]

    polPower :: [Double] -> Int -> [Double]
    polPower f 0 = [1]
    polPower f n = multPoly f (polPower f (n-1))

    polDivCoeff :: [Double] -> [Double] -> Double
    polDivCoeff f g = head f / head g

    polDivRemain :: [Double] -> [Double] -> Double -> [Double]
    polDivRemain f g coeff = -- g - f
        if length f < length g then f
        else
        sumVectors (map (*(-coeff)) g) f

    f :: [Double]
    f = [2, 3, 1]
    g :: [Double]
    g = [1, -1]

    polDiv :: [Double] -> [Double] -> ([Double], [Double], [Double])
    -- Dividing two polynomials, f, g. Returns division + reminader (num + denom)
    -- Note: MUST send in division as flipped polynomial ([1, 2, 3] => x^2 + 2x + 3)
    polDiv f g = -- f / g
        if length f < length g then ([], f, g)
        else
            let coeff = polDivCoeff f g
                remain = polDivRemain f g coeff
                (x, y, z) = polDiv (tail remain) g
            in (coeff:x, y, g)


    -- computes pol rep for f (g (x)), given approximations f and g
    fnComposition :: [Double] -> [Double] -> [Double]
    fnComposition f g =
        let zipped = zip f [0..(length f)]
            composed = map (\ (x, y) -> map (*x) (polPower g y))  zipped -- x is coefficient. y is order of polynomial. we need to multiply x by exp of g
        in foldl (\x y -> sumVectors x y) [] composed
