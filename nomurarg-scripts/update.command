#!/bin/bash
cd -- "$(dirname "$BASH_SOURCE")"

git pull nomura master
git pull origin master