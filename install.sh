#!/bin/bash
virtualenv -p$(which python3) env || exit
source env/bin/activate || exit
git clone https://github.com/mtholder/peyotl.git || exit

cd peyotl || exit
git checkout -b taxalotl origin/taxalotl || exit
pip install -r requirements.txt || exit
python setup.py develop || exit
cd - 2>/dev/null || exit

python setup.py develop || exit

echo
echo 'Taxalotl installed, but you will still need to run:'
echo '    source env/bin/activate'
echo 'in your shell before you can access it.'

usershome=$(readlink -f ~)
fptd="${usershome}/.opentreeoflife/taxalotl/"
if ! test -d "${fptd}"
then
    mkdir -p "${fptd}" || exit
fi

if ! test -f "${fptd}/taxalotl.conf"
then
    cat taxalotl.conf.example | sed -e "s:ABS-PATH-TO-YOUR-HOME-DIR:${usershome}:" > "${fptd}/taxalotl.conf" || exit
fi

echo
echo "You almost certainly want to run something like:"
echo "   export PEYOTL_CONFIG_FILE=\"$PWD/recommended-peyotl-conf.ini\""
echo "in your bash shell (or the equivalent if you do not use bash for some reason)"
echo "so that you can see the logging messages."

