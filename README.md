# MaNGA-BCGs
This pipeline is designed for removing MaNGA systematics on MPL4/5.
## Get started
The script requires [sep](https://sep.readthedocs.io/en/v1.0.x/) and [swarp](http://www.astromatic.net/software/swarp).
  
## Usage

`bcg_selection.py` is used to select BCGs and write them in the correct input format to feed into further scripts.

`sdss.py` extracts satellites and downloads jpegs of the selected BCGs.

`bcg_xyvviar.py` is used to stack all the BCG velocity maps and masks.

`check_bcg_xyvviar_mpl5.py` checks individual BCGs.

