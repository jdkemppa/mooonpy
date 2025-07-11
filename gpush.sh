# This script is meant for quick Git pushes for Josh's development purposes
if [ $# -eq 0 ]; then
  # Do not use this option if you are worried others have access to your machine!
  echo "No arguments provided. Using internally set message"
  MESSAGE="gpush.sh quick - JDK"
else
  MESSAGE=$1
  echo "Setting message as '$MESSAGE'"
fi

cd docs/
make clean
cd ../

git pull
git add .
git commit -m "$MESSAGE"
git push origin main
