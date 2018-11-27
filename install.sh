#!/bin/bash
virtualenv venv || exit
source venv/bin/activate || exit
git clone https://github.com/mtholder/peyotl.git || exit

cd peyotl || exit
git checkout taxalotl || exit
pip install -r requirements.txt || exit
python setup.py develop || exit
cd - 2>/dev/null || exit

python setup.py develop || exit

echo
echo 'Taxalotl installed, but you will still need to run:'
echo '    source venv/bin/activate'
echo 'in your shell before you can access it.'

if ! test -f taxalotl.conf
then
    cp taxalotl.conf.example taxalotl.conf || exit
    echo 'You will also need to edit "taxalotl.conf" to choose locations for the data'
fi

echo
echo "You almost certainly want to run something like:"
echo "   export PEYOTL_CONFIG_FILE=\"$PWD/recommended-peyotl-conf.ini\""
echo "in your bash shell (or the equivalent if you do not use bash for some reason"
echo "so that you can see the logging messages."

