

# This is a script used to setup a virtual environment that includes the
# mevlib package in path and enter that virtual environment to develop the
# package.  (The package will not be cached in the way most packages are.)


# Because we want this to be executed in the current shell rather than
# another bash process, you must 'source' this file.  In other words,
#   $ bash ./start_developing.sh
# will not work, but
#   $ source ./start_developing.sh
# will.


python3 -m venv virtualenv
source ./virtualenv/bin/activate
pip3 install wheel      # not doing this first can cause pip errors
pip3 install -r requirements.txt
python3 setup.py develop


echo "You are now ready to develop the package and/or run scripts in this    "
echo "package without adjusting the path.  Do leave this virtual environment,"
echo "use 'deactivate'."


# the following command will leave the virtual environment
#deactivate


