---
layout: post
title:  "Choose Colors for Scientific Data Visualization"
date:   2021-04-06 20:35:58 +0800
usemathjax: true
categories: jekyll update
---


## Some general suggestion
- [Choosing color palettes for scientific figures](https://onlinelibrary.wiley.com/doi/full/10.1002/rth2.12308)
- [Color in Scientific Figures](https://www.aje.com/arc/Using-Color-in-Figures/)
- [Design Tips for Scientists](https://cns.utexas.edu/images/CNS/Deans_Office/Communications/Files/design-tips-for-scientists_GUIDE.pdf)
- 2014, Plos Computational Biology, [Ten Simple Rules for Better Figures](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833)


## Background in color space
1. RGB color space 
   - Mixture of red, green and blue color

2. HSB color space
   - **Hue** (色相)
   - **Saturation** (饱和度)
   - **Brightness** (亮度)， some tines also called **value** or **lightness**

- See <https://www.codeproject.com/Articles/1202772/Color-Topics-for-Programmers>
- python for color space mapping
- For detailed implementation of such conversion see:
  - <https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both>

```python
import colorsys
# convert rgb to hsv 
colorsys.rgb_to_hsv(0.2, 0.4, 0.4)
# (0.5, 0.5, 0.4)
colorsys.hsv_to_rgb(0.5, 0.5, 0.4)
# (0.2, 0.4, 0.4)

```

3. [HSLuv](https://www.hsluv.org/) 
- A alternative to `HSL`, claimed to be more human friendly

## Some useful tricks
- Given rgb colors:
  1. Change the lightness / saturation 
    - The easist way may be convert `rgb` values to `hsv` values, change `s`/`v`, than map back to rgb value
  2. Take complementary color 
    - Convert to `hsv`, then set `h=(h+0.5)%1`
    - Also see <https://stackoverflow.com/questions/40233986/python-is-there-a-function-or-formula-to-find-the-complementary-colour-of-a-rgb> 
  3. Interpolate from one color to another
    - We shall still convert `rgb` to `hsv`, then `h` ...
    - See discussion in <https://stackoverflow.com/questions/13488957/interpolate-from-one-color-to-another>

## Color map for visualization
### Color Map Requirements
- According to [Diverging Color Maps for Scientific Visualization](http://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf)
  – The map yields images that are aesthetically pleasing.
  – The map has a maximal perceptual resolution.
  – Interference with the shading of 3D surfaces is minimal.
  – The map is not sensitive to vision deficiencies.
  – The order of the colors should be intuitively the same for all people.
  – The perceptual interpolation matches the underlying scalars of the map
### Claasification of color maps
- qualitative / nominal color maps:  discrete, unordered classes
-  sequential / ordinal / saturation color maps: hue is nearly fixed, (nearly monochromatic), difference in saturation and lightness indicate numericial difference
- diverging / ratio / bipolar / doubleended  color maps:  two major color components
- Cyclic: change in lightness of two different colors that meet in the middle and beginning/end at an unsaturated color

### Collection of perceptually accurate colormaps
  - <https://colorcet.holoviz.org/>
  - [colorcet](https://colorcet.com/), <https://github.com/holoviz/colorcet>
  - <https://www.kennethmoreland.com/color-maps/>
  - <https://github.com/1313e/CMasher>

### Caveats
- <https://arxiv.org/abs/1509.03700>
- **Don't** use color map like `jet` or `hot`
  - 2020, Nature Communication, [The misuse of colour in science communication](https://www.nature.com/articles/s41467-020-19160-7)

## Examples

### python
- Spectral, coolwarm, RdBu, vlag


## Resources
### Find some palette 
  - <https://pattern-library.economist.com/color.html> 
### R packages
- [ggsci](https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html)
- [Wes Anderson Palettes](https://github.com/karthik/wesanderson)
- [GeoDataViz-Toolkit](https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/tree/master/Colours)



