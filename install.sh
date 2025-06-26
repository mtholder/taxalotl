#!/bin/bash
set -x

cd .. || exit
python3 -mvenv env || exit
source env/bin/activate || exit

if ! test -d peyutil ; then
    git clone git@github.com:OpenTreeOfLife/peyutil.git || exit
fi
cd peyutil || exit
git checkout -b modern-taxalotl origin/modern-taxalotl || exit
pip install -r requirements.txt || exit
pip install -e . || exit
cd - 2>/dev/null || exit

if ! test -d peyotl ; then
    git clone git@github.com:OpenTreeOfLife/peyotl.git || exit
fi
cd peyotl || exit
git checkout -b modern-taxalotl origin/modern-taxalotl || exit
pip install -r requirements.txt || exit
pip install -e . || exit
cd - 2>/dev/null || exit

cd taxalotl || exit
pip install -r requirements.txt || exit
pip install -e . || exit

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

