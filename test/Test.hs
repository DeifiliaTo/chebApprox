module Main (main) where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable, (*~), Time )
import Diagrams.Backend.CmdLine

import System.IO.Unsafe

import ChebyshevApproximations

denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
 mainRender ( DiagramOpts (Just 900) (Just 700) fn
            , DiagramLoopOpts False Nothing 0
            )

chart :: String ->
        [[(Double, Double)]] ->
        Renderable ()
chart t obss = toRenderable layout
 where

   actual :: [(Double, Double)] -> Colour Double -> String -> PlotLines Double Double
   actual x c t = plot_lines_values .~ [x]
                $ plot_lines_style  . line_color .~ opaque c
                $ plot_lines_style  . line_width .~ 1.0
                $ plot_lines_title .~ t
                $ def

   cs = cycle [blue, green, red]
   ts = map show [1..]

   actuals' :: [PlotLines Double Double]
   actuals' = zipWith3 actual obss cs ts

   layout = layout_title .~ t
          $ layout_plots .~ (map toPlot actuals')
          $ layout_y_axis . laxis_title .~ "Value"
          $ layout_y_axis . laxis_override .~ axisGridHide
          $ layout_x_axis . laxis_title .~ "Time"
          $ layout_x_axis . laxis_override .~ axisGridHide
          $ def

diag :: String ->
       [[(Double, Double)]] ->
       Diagram Cairo
diag t xss =
 fst $ runBackend denv (render (chart t xss) (600, 500))

main :: IO ()
main = return ()
