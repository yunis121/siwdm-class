set -e

echo "This is the installation routine for the siwdm-class module. It will unzip CLASS v2.7.2, replace the modified files and compile the source code. For any installation problems we refer to the handy CLASS intallation guide at https://github.com/lesgourg/class_public"

echo "Unzipping"

unzip class_public-2.7.zip

echo "Copying source"

cp -R ./include ./class_public-2.7

cp -R ./source ./class_public-2.7

cp explanatory_plus_siwdm.ini ./class_public-2.7

cp explanatory_plus_siwdm.pre ./class_public-2.7

cp explanatory_plus_siwdm_reltime.dat ./class_public-2.7

cp explanatory_plus_siwdm_reltime_params.txt ./class_public-2.7

echo "Rebuilding CLASS..."

cd ./class_public-2.7

make clean

make

echo "Rebuild Complete!"

echo "SIWDM module sucessfully installed on top of CLASS v2.7.2 at"$PWD

echo "Open explanatory_plus_siwdm.ini to get started and run the program using ./class explanatory_plus_siwdm.ini explanatory_plus_siwdm.pre "
