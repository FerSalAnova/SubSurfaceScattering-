# Subsurface Scattering in Nori (with KD-Tree Acceleration)

This repository implements **subsurface scattering (SSS)** for the [Nori](https://github.com/wjakob/nori) physically-based renderer using **point-based preprocessing** and a **KD-tree search** to significantly improve execution time.

Our method reduces computational complexity to approximately **O(n³)** by organizing sample points efficiently in a spatial data structure.

---

## Features

- Implements **point-based subsurface scattering**
- Uses a **KD-tree** for efficient lookup and evaluation
- Achieves improved performance over naive approaches
- Customizable via Nori XML scene descriptions

---

## 🛠Installation & Setup

1. **Install Nori**  
   Follow the official instructions from the [Nori repository](https://github.com/wjakob/nori).

2. **Clone or copy this repository** into your local Nori directory.

3. **Merge the following folders** into the corresponding Nori directories:
/include/ → Add to Nori's include path
/src/ → Add to Nori's src directory

4. **Rebuild Nori** after merging:

```bash
mkdir build
cd build
cmake ..
make -j
```
## Resulting Image

Here is an example of the rendered scene with subsurface scattering enabled:

![SSS Render](Figures/Final_scene.png)
