#! /bin/sh

echo VERSION=\"`git describe`\" > version.py
echo TIMESTAMP=\"`git log -1 --date=iso --pretty="format:%cd"`\" >> version.py
