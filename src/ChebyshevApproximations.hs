module ChebyshevApproximations
where
    import Prelude
    
  
    -- Returns value of chebyshev zero
    computeChebNode :: (Floating a, Integral b) => b -> b -> a
    computeChebNode n k = 
        cos ((2*fromIntegral(n)+1-2*fromIntegral(k))*pi/(2*fromIntegral(n)+2))

    -- Creates list of chebyshev zeros
    chebNodes :: (Floating a, Integral b) => b -> [a]
    chebNodes n =
        let nodes = [0..n] in
            map (\x -> computeChebNode n x) (nodes)

    -- Takes in an order. Returns list of chebyshev polynomials
    -- Takes in an order. Returns list of chebyshev polynomials
    chebPol :: (Floating a, Integral b) => b -> [[a]]
    chebPol 0 = [[1.0]]
    chebPol 1 = [[0.0, 1.0], [1.0]]
    chebPol n =
        let prevResult = chebPol (n-1) in
        let multTwo = map (\x -> map (*2.0) x) prevResult in
            let firstTerm = map (\x -> 0:x) multTwo in
                let subtractTerm = sumVectors (head (firstTerm)) ((map (*(-1)) (head (tail prevResult)))) in
                    subtractTerm:(prevResult)


    -- calcCoeffs :: ([Double] -> [Double]) -> Double -> Double
    -- will be much easier w/ vector operation
    -- Adds up all the elements in a vector
    -- TODO: Replace with vector ops/accelerate library
    sumList :: (Floating a) => [a] -> a
    sumList lst = foldl (+) 0.0 lst

    -- Computes f (x_k) * cos (x_k)
    computeProduct :: (Floating a, Integral b) => (a -> a) -> b -> b -> b -> a
    computeProduct f n j k =
        f ((2*fromIntegral(k)+1)/(2*fromIntegral(n)+2)*pi) * cos (fromIntegral (j)*(2*fromIntegral(k)+1)/(2*fromIntegral(n)+2)*pi)

    -- Calculates c_j (coefficient)
    chebCoeff :: (Floating a, Integral b) => (a -> a) -> b  -> b -> [a]
    chebCoeff f n iter =
        let nodes = chebNodes n in
            if iter == 0 then [1/(1*(fromIntegral(n)+1))*sumList (map f nodes)]
            else
                2.0/(fromIntegral(n)+1.0)*sumList(
                    map (\x ->
                        --computeProduct f n iter x
                        f (computeChebNode n x)*cos (fromIntegral(iter) *(2*fromIntegral(n)+1-2*fromIntegral(x))/(2*fromIntegral(n)+2)*pi)
                    ) [0..n]
                ):(chebCoeff f n (iter-1))

    -- function to add two polynomials together
    sumVectors :: (Floating a) => [a] -> [a] -> [a]
    sumVectors p1 p2 =
        if (length p1 >= length p2)
        then zipWith (+) p1 (p2 ++ repeat 0)
        else sumVectors p2 p1


    -- takes in list of coefficients, and cheb polynomials, outputs polynomial approximation
    -- Used for c_j * T_j


    chebf :: (Floating a, Integral b) => (a -> a) -> b -> [a]
    chebf f n =
        let coeffs = chebCoeff f n n
            chebPols = chebPol n
            zipped = zip coeffs chebPols
            mapped = map (\(x, y) -> map (*x) y) zipped
            in foldl (\x y -> sumVectors x y) [] mapped

    polCalc :: (Floating a, Eq a, Integral b) => [a] -> a -> b -> a
    polCalc (c:coeffs) x order =
        if coeffs == [] then
            c * x^order
        else
            c * x^order + polCalc coeffs x (order + 1)



    -- Example functions


    h :: Double -> Double
    h x = sin x + exp (cos x) - x*x*x


    -- Newtonian approach
    findDiff :: (Floating a, Integral b) => a -> a -> b -> a
    findDiff x y  order = (y-x)/fromIntegral(order)

    -- Given [1, 2, 4, 8], order 1 -> [1, 2, 4]
    divDiff :: (Floating a, Eq a, Integral b) =>  [a] -> b -> [a]
    divDiff (x:xs) order =
        if xs == [] then []
        else
            let b = head xs in
                (findDiff x b order):(divDiff xs order)

    -- Given [1, 2, 4, 8] --> [[1, 2, 4, 8], [1, 2, 4], [0.5, 1], [1/6]]
    divDiffList :: (Floating a, Eq a, Integral b) => [a] -> b -> [[a]]
    divDiffList (x:xs) order =
       if length (x:xs) == 2 then [divDiff (x:xs) (order)]
       else
            let lst = (divDiff (x:xs) order) in
                lst:(divDiffList lst (order+1))

    newtonCoeffs :: (Floating a, Eq a) => [a] -> [a]
    newtonCoeffs lst =
        (head lst):(map head (divDiffList lst 1))

    multiplyByX :: (Floating a) => [a] -> [a]
    multiplyByX p = 0:p

    multPoly :: (Floating a) => [a] -> [a] -> [a]
    multPoly [] p2 = []
    multPoly (p:p1) p2 =
        let pTimesP2 = map (*p) p2 in
            let xTimesP1Timesp2 = multiplyByX (multPoly p1 p2) in
                sumVectors pTimesP2 xTimesP1Timesp2

    -- Generate list of polynomials to be multiplied (Newtonian form)
    newtonPolMult :: (Floating a) => [a] -> [[a]]
    newtonPolMult nodes =
        map (\x -> [-1*x, 1]) nodes

    polMultFoldable :: (Floating a) => [[a]] -> [[a]]
    polMultFoldable lst =
        [1]:scanl1 (\x y -> multPoly x y) lst

    fn :: Double -> Double
    fn x = 2**x

    newtonApprox :: (Floating a, Integral b, Eq a, Enum a) => (a -> a) -> b -> [a]
    newtonApprox f n =
        --let newtonNodes = chebNodes n
        let newtonNodes = [0.5, 0.6..2]
            coeffs = newtonCoeffs (map f newtonNodes)
            pairedList = zip coeffs (polMultFoldable (newtonPolMult newtonNodes))
            mapped = map (\(x, y) -> map (*x) y) pairedList
            in foldl (\x y -> sumVectors x y) [] mapped

    printErrLag :: (Floating a, Integral b, Eq a, Enum a) => (a -> a) -> b -> [a]
    printErrLag f order =
        let soln = chebf f order in
            map (\x -> (polCalc soln x 0) - f x ) [-1, -0.99..1]

    printErrNew :: (Floating a, Integral b, Eq a, Enum a) => (a -> a) -> b -> [a]
    printErrNew f order =
        let soln = newtonApprox f order in
            map (\x -> (polCalc soln x 0) - f x ) [0.5, 1..2]

    polPower :: (Floating a, Integral b) => [a] -> b -> [a]
    polPower f 0 = [1]
    polPower f n = multPoly f (polPower f (n-1))

    polDivCoeff :: (Floating a) => [a] -> [a] -> a
    polDivCoeff f g = head f / head g

    polDivRemain :: (Floating a) => [a] -> [a] -> a -> [a]
    polDivRemain f g coeff = -- g - f
        if length f < length g then f
        else
        sumVectors (map (*(-coeff)) g) f

    f :: [Double]
    f = [2, 3, 1]
    g :: [Double]  
    g = [1, -1]

    polDiv :: (Floating a) => [a] -> [a] -> ([a], [a], [a])
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
    fnComposition :: (Floating a) => [a] -> [a] -> [a]
    fnComposition f g =
        let zipped = zip f [0..(length f)]
            composed = map (\ (x, y) -> map (*x) (polPower g y))  zipped -- x is coefficient. y is order of polynomial. we need to multiply x by exp of g
        in foldl (\x y -> sumVectors x y) [] composed

    envelope :: (Floating a, Eq a, Ord a) => [a] -> [Int] -> [(a, Int)]
    envelope (c:coeffs) (i:ndices) = 
        if (coeffs) == [] then [(c, i)]
        else 
            let env = envelope coeffs ndices in
                if c > (fst (head env)) then ((maximum (c:coeffs)), i):(env)
                else env

    
    tol :: (Floating a) => a
    tol = 1e-15

    toInt :: (Integral a, RealFrac a) => a -> Int
    toInt x = round (x)

    --plateau :: (Floating a, Ord a, Integral b, Fractional b) => [a] -> b -> b -> b
    plateau :: (Floating a, Ord a) => [(a, Int)] -> Int -> Int -> Int
    plateau env j n =
        let j2 = round ((1.25*fromIntegral(j))+5.0) in
            if j2 <= (length (env)-1) then
                if (fst(env!!j2)/fst(env!!j)) >=  (3.0*(1.0-log(fst (env!!j))/log(tol)))
                    then j-1
                else
                    plateau env (j+1) n
            else snd (env!!j)

    lst :: [Double]
    lst = [1, 4.0, 3.0]

    extractEnv :: (Floating a) => [(a, Int)] -> [a]
    extractEnv lst = map (\x -> fst x) lst
    -- findIndex :: (Floating a, Eq a) => [a] -> Int
    --findIndex coeffs =
    --  toIntplateau (envelope coeffs [0..(length coeffs -1)]) 0 (length coeffs)