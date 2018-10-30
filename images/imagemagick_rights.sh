#!/bin/bash
# ImageMagick has very strict read policies since version 6 (6?) and disallows reading PDFs.
# To be able to convert PDFs to PNGs, we have to change these policy settings.
POLICY_PATH="/etc/ImageMagick-6/policy.xml"
perl -pe 's/(?<=policy domain="coder" rights=")none(?=" pattern="PDF")/read/g' $POLICY_PATH > policy.xml
sudo mv policy.xml /etc/ImageMagick-6/policy.xml
