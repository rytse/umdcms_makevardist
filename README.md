# umdcms_makevardist

Scripts for looking for interesting cuts of the CMS data.

## Setup

In order to run this on the UMD T3 cluster you must do the following:
* Set up your CMSSW environment
* Set up the UMD CMS group's analysis library `WG_Analysis`
* Install missing python dependencies (locally)

The `WG_Analysis` repo kindly provides a script that sets up the correct CMSSW environment, so all we need to do is to be sure to pull this repo recusively, run that script, and then install our python dependencies.

`$ git clone --recursive git@github.com:rytse/umdcms_makevardist.git` if you prefer SSH, or

`$ git clone --recursive https://github.com/rytse/umdcms_makevardist.git` if you prefer HTTPS

Then run the appropriate setup script from `WG_Analysis` (`setup.sh` for bash, `setup.csh` for heathens, etc.)

Following this, use pip to install `uncertainties` to your local home directory
`$ pip install --user uncertainties`

and be sure to add the appropriate directory to your python path.
