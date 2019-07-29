module Main (main) where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable, (*~), Time )
import Diagrams.Backend.CmdLine
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.LLVM.Native as CPU
import qualified Prelude as P
import Test.HUnit


import System.IO.Unsafe

import ChebApproxAcc
import ChebTypingAcc


denv :: DEnv A.Double
denv = unsafePerformIO P.$ defaultEnv vectorAlignmentFns 600 500

displayHeader :: P.FilePath -> Diagram B -> P.IO ()
displayHeader fn =
 mainRender ( DiagramOpts (P.Just 900) (P.Just 700) fn
            , DiagramLoopOpts P.False P.Nothing 0
            )

chart :: P.String ->
        [[(A.Double, A.Double)]] ->
        Renderable ()
chart t obss = toRenderable layout
 where

   actual :: [(A.Double, A.Double)] -> Colour A.Double -> P.String -> PlotLines A.Double A.Double
   actual x c t = plot_lines_values .~ [x]
                P.$ plot_lines_style  P.. line_color .~ opaque c
                P.$ plot_lines_style  P.. line_width .~ 1.0
                P.$ plot_lines_title .~ t
                P.$ def

   cs = P.cycle [blue, green, red]
   ts = P.map P.show [1..]

   actuals' :: [PlotLines A.Double A.Double]
   actuals' = P.zipWith3 actual obss cs ts

   layout = layout_title .~ t
          P.$ layout_plots .~ (P.map toPlot actuals')
          P.$ layout_y_axis P.. laxis_title .~ "Value"
          P.$ layout_y_axis P.. laxis_override .~ axisGridHide
          P.$ layout_x_axis P.. laxis_title .~ "Time"
          P.$ layout_x_axis P.. laxis_override .~ axisGridHide
          P.$ def

diag :: P.String ->
       [[(A.Double, A.Double)]] ->
       Diagram Cairo
diag t xss =
 P.fst P.$ runBackend denv (render (chart t xss) (600, 500))

main = 
       let 
              f = ChebApproxAcc.f
              res = ChebApproxAcc.chebf f 8
              residual = A.toList (CPU.run P.$ ChebApproxAcc.approxError f res)
              {- x = ChebTyping.x 
              f = fromCheb (ChebTyping.cos x)
              cheb = map (\x -> (polCalc f x 0)) [-1, -0.99..1]
              residual = map (\x -> (polCalc f x 0) Prelude.- (Prelude.cos x)) [-1, -0.99..1] 
              cheb20 = chebf Prelude.cos 20
              cheb20plot = map (\x -> (polCalc (cheb20) x 0) Prelude.- Prelude.cos x ) [-1, -0.99..1] -}
       in 
       do
              --renderableToFile def "cos-cheb" (chart "cos-cheb" ([zip [-1, -0.99..1] cheb]))
              --renderableToFile def "cos-residual" (chart "cos-residual" ([P.zip [-1, -0.9..1] residual]))
              --renderableToFile def "cos-residual-n=20" (chart "cos-residual-n=20" ([zip [-1, -0.99..1] cheb20plot]))