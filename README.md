# MaNGA-BCGs
This pipeline is designed for removing MaNGA systematics on MPL4/5.
## Getting started
The script requires [sep](https://sep.readthedocs.io/en/v1.0.x/) and [swarp](http://www.astromatic.net/software/swarp).
  
## Usage

`bcg_selection.py` is used to select BCGs and write them in the correct input format to feed into further scripts.

`sdss.py` extracts satellites and downloads jpegs of the selected BCGs.

`bcg_xyvviar.py` is used to stack all the BCG velocity maps and masks.

<img src="https://cloud.githubusercontent.com/assets/6698757/19276763/7361e79e-8fa5-11e6-87fc-58792aee83da.png" width="1000">

`check_bcg_xyvviar_mpl5.py` checks individual BCGs. Use `make` to plot all individual BCGs.

<img src="https://cloud.githubusercontent.com/assets/6698757/19276683/2abc8544-8fa5-11e6-9204-f8c0f508db04.png" width="1000">
