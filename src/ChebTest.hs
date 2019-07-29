{-# L GUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DeriveAnyClass  #-}
{-# LANGUAGE DeriveGeneric   #-}
{-# LANGUAGE PatternSynonyms #-}

module  ChebTest
where
    import Prelude as P
    import Data.Array.Accelerate as A
    import Data.Array.Accelerate.Debug as A
    import Data.Array.Accelerate.Interpreter as I
    import Data.Array.Accelerate.LLVM.Native as CPU
    import ChebTypingAcc
    import ChebApproxAcc
    import ChebMath
    import Test.HUnit

    -- [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    arr1 :: Acc (Vector Double)
    arr1 = use (fromList (Z :. 10) [0..])
    
    arr1' :: Vector Double
    arr1' = CPU.run $ arr1

    -- [1, 3, 2]
    arr2 :: Acc (Vector Double)
    arr2 = use (fromList (Z :. 3) [1, 3, 2])

    arr2' :: Vector Double
    arr2' = CPU.run $ arr2

    -- [1]
    arr3 :: Acc (Vector Double)
    arr3 = use (fromList (Z :. 1) [1])
    
    arr3' :: Vector Double
    arr3' = CPU.run $ arr3

    -- [0, 0, -1]
    arr4 :: Acc (Vector Double)
    arr4 = use (fromList (Z :. 3) [0, 0, -1])

    arr4' :: Vector Double
    arr4' = CPU.run $ arr4


    addTest0 = TestCase 
        (
            assertEqual "[0..9] + [1, 3, 2]"
            (fromList (Z :. 10) [1, 4, 4, 3, 4, 5, 6, 7, 8, 9]) 
            (CPU.run $ sumVectors arr1 arr2)
        )
    
    addTest1 = TestCase 
        (
            assertEqual "[1, 3, 2] + [0..9]" 
            (fromList (Z :. 10) [1, 4, 4, 3, 4, 5, 6, 7, 8, 9])
            (CPU.run $ sumVectors arr2 arr1)
        )

    addTest2 = TestCase 
        (
            assertEqual "[0..9] + [1]" 
            (fromList (Z :. 10) [1, 1, 2, 3, 4, 5, 6, 7, 8, 9])
            (CPU.run $ sumVectors arr1 arr3)
        )

    addTest3 = TestCase 
        (
            assertEqual "[1] + [0..9]" 
            (fromList (Z :. 10) [1, 1, 2, 3, 4, 5, 6, 7, 8, 9])
            (CPU.run $ sumVectors arr3 arr1)
        )

    addTest4 = TestCase 
        (
            assertEqual "[0..9] + [0, 0, -1]"
            (fromList (Z :. 10) [0, 1, 1, 3, 4, 5, 6, 7, 8, 9]) 
            (CPU.run $ sumVectors arr1 arr4)
        )

    addTest5 = TestCase 
        (
            assertEqual "[0, 0, -1] + [0..9]" 
            (fromList (Z :. 10) [0, 1, 1, 3, 4, 5, 6, 7, 8, 9])
            (CPU.run $ sumVectors arr4 arr1)            
        )
    
    addTest6 = TestCase 
        (
            assertEqual "[1, 3, 2] + [1]" 
            (fromList (Z :. 3) [2, 3, 2])
            (CPU.run $ sumVectors arr2 arr3)
        )

    addTest7 = TestCase 
        (
            assertEqual "[1] + [1, 3, 2]" 
            (fromList (Z :. 3) [2, 3, 2])
            (CPU.run $ sumVectors arr3 arr2)
        )

    addTest8 = TestCase 
        (
            assertEqual "[1, 3, 2] + [0, 0, -1]" 
            (fromList (Z :. 3) [1, 3, 1])
            (CPU.run $ sumVectors arr2 arr4)
        )

    addTest9 = TestCase 
        (
            assertEqual "[0, 0, -1] + [1, 3, 2]" 
            (fromList (Z :. 3) [1, 3, 1])
            (CPU.run $ sumVectors arr4 arr2)
        )
    
    addTest10 = TestCase 
        (
            assertEqual "[1] + [0, 0, -1]" 
            (fromList (Z :. 3) [1, 0, -1])
            (CPU.run $ sumVectors arr3 arr4)
        )

    addTest11 = TestCase 
        (
            assertEqual "[0, 0, -1] + [1]" 
            (fromList (Z :. 3) [1, 0, -1])
            (CPU.run $ sumVectors arr4 arr3)
        )
    
    addTests = TestList 
        [
            TestLabel "add test 0" addTest0,
            TestLabel "add test 1" addTest1,
            TestLabel "add test 2" addTest2,
            TestLabel "add test 3" addTest3,
            TestLabel "add test 4" addTest4,
            TestLabel "add test 5" addTest5,
            TestLabel "add test 6" addTest6,
            TestLabel "add test 7" addTest7,
            TestLabel "add test 8" addTest8,
            TestLabel "add test 9" addTest9,
            TestLabel "add test 10" addTest10,
            TestLabel "add test 11" addTest11
        ]

    multTest0 = TestCase 
        (
            assertEqual "[0..9] * [1, 3, 2]"
            (fromList (Z :. 12) [0, 1, 5, 11, 17, 23, 29, 35, 41, 47, 43, 18])
            (CPU.run $ multPoly arr1 arr2)
        )
    
    multTest1 = TestCase 
        (
            assertEqual "[1, 3, 2] * [0..9]" 
            (fromList (Z :. 12) [0, 1, 5, 11, 17, 23, 29, 35, 41, 47, 43, 18])
            (CPU.run $ multPoly arr2 arr1)
        )

    multTest2 = TestCase 
        (
            assertEqual "[0..9] * [1]" 
            (fromList (Z :. 10) [0..])
            (CPU.run $ multPoly arr1 arr3)
        )

    multTest3 = TestCase 
        (
            assertEqual "[1] * [0..9]" 
            (fromList (Z :. 10) [0..])
            (CPU.run $ multPoly arr3 arr1)
        )

    multTest4 = TestCase 
        (
            assertEqual "[0..9] * [0, 0, -1]"
            (fromList (Z :. 12) [0, 0, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9]) 
            (CPU.run $ multPoly arr1 arr4)
        )

    multTest5 = TestCase 
        (
            assertEqual "[0, 0, -1] * [0..9]" 
            (fromList (Z :. 12) [0, 0, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9]) 
            (CPU.run $ multPoly arr4 arr1)            
        )
    
    multTest6 = TestCase 
        (
            assertEqual "[1, 3, 2] * [1]" 
            (fromList (Z :. 3) [1, 3, 2])
            (CPU.run $ multPoly arr2 arr3)
        )

    multTest7 = TestCase 
        (
            assertEqual "[1] * [1, 3, 2]" 
            (fromList (Z :. 3) [1, 3, 2])
            (CPU.run $ multPoly arr3 arr2)
        )

    multTest8 = TestCase 
        (
            assertEqual "[1, 3, 2] * [0, 0, -1]" 
            (fromList (Z :. 5) [0, 0, -1, -3, -2])
            (CPU.run $ multPoly arr2 arr4)
        )

    multTest9 = TestCase 
        (
            assertEqual "[0, 0, -1] * [1, 3, 2]" 
            (fromList (Z :. 5) [0, 0, -1, -3, -2])
            (CPU.run $ multPoly arr4 arr2)
        )
    
    multTest10 = TestCase 
        (
            assertEqual "[1] * [0, 0, -1]" 
            (fromList (Z :. 3) [0, 0, -1])
            (CPU.run $ multPoly arr3 arr4)
        )

    multTest11 = TestCase 
        (
            assertEqual "[0, 0, -1] * [1]" 
            (fromList (Z :. 3) [0, 0, -1])
            (CPU.run $ multPoly arr4 arr3)
        )
    
    multTests = TestList 
        [
            TestLabel "mult test 0" multTest0,
            TestLabel "mult test 1" multTest1,
            TestLabel "mult test 2" multTest2,
            TestLabel "mult test 3" multTest3,
            TestLabel "mult test 4" multTest4,
            TestLabel "mult test 5" multTest5,
            TestLabel "mult test 6" multTest6,
            TestLabel "mult test 7" multTest7,
            TestLabel "mult test 8" multTest8,
            TestLabel "mult test 9" multTest9,
            TestLabel "mult test 10" multTest10,
            TestLabel "mult test 11" multTest11
        ]
    
    scalarMultTest0 = TestCase
        (
            assertEqual "0 * [0..n]"
            (CPU.run $ ( A.fill (lift (Z:. (10::Exp Int))) 0))
            (CPU.run $ fromCheb (0 *^ (toCheb arr1)))
        )
    
    scalarMultTest1 = TestCase
        (
            assertEqual "[0..n] * 0"
            (CPU.run $ ( A.fill (lift (Z:. (10::Exp Int))) 0))
            (CPU.run $ fromCheb ((toCheb arr1) ^* 0))
        )

    scalarMultTest2 = TestCase
        (
            assertEqual "1 * [1]"
            (fromList (Z :. 1) [1])
            (CPU.run $ fromCheb (1 *^ (toCheb arr3)))
        )
    
    scalarMultTest3 = TestCase
        (
            assertEqual "[1] * 1"
            (fromList (Z :. 1) [1])
            (CPU.run $ fromCheb ((toCheb arr3) ^* 1))
        )
    
    scalarMultTests = TestList
        [
            TestLabel "scalar mult test 0" scalarMultTest0,
            TestLabel "scalar mult test 1" scalarMultTest1,  
            TestLabel "scalar mult test 2" scalarMultTest2,
            TestLabel "scalar mult test 3" scalarMultTest3 
        ]