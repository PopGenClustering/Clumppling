python -m clumppling \
-i input/diffInd/ind200 \
-o output/diffInd/ind200 \
-f structure --cd_default 0 
python -m clumppling \
-i input/diffInd/ind100 \
-o output/diffInd/ind100 \
-f structure --cd_default 0 
python -m clumppling \
-i input/diffInd/ind50 \
-o output/diffInd/ind50 \
-f structure --cd_default 0 
python -m clumppling.diffInd \
-i output/diffInd \
-o output/diffInd 