#!/bin/bash
# This script lets Travis CI deploy the FORD generated documentation website

set -ev

if [ ! "$TRAVIS" ]; then
    echo "Documentation can only be deployed by Travis CI"
    exit 0
fi

if [ "$TRAVIS_BRANCH" != "master" ]; then
  echo "Skipping documentation deployment"
  echo "Only applicable to master branch"
  exit 0
fi

if [ "$TRAVIS_PYTHON_VERSION" != "3.4" ]; then
  echo "Skipping documentation deployment"
  echo "Only use Python 3.4"
  exit 0
fi


if [ "$GCC_VERSION" != "7" ]; then
  echo "Skipping documentation deployment"
  echo "Only use gcc7"
  exit 0
fi

if [ "$BUILD_TYPE" != "Release" ]; then
  echo "Skipping documentation deployment"
  echo "Only use Release build"
  exit 0
fi

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    echo "Skipping documentation deployment"
    echo "Not applicable to pull requests"
    exit 0
fi

REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
echo $REPO
echo $SSH_REPO

eval `ssh-agent -s`
ssh-add .deploy_key

git clone --branch=gh-pages $REPO gh-pages

cd gh-pages
rm -rf css favicon.png fonts index.html interface \
   js lists media module page proc program search.html \
   sourcefile src tipuesearch type
cp -r "$TRAVIS_BUILD_DIR"/docs/html/* .
git add *
git commit -am "Documentation updated by travis job $TRAVIS_JOB_NUMBER for commits $TRAVIS_COMMIT_RANGE"
git push $SSH_REPO gh-pages
