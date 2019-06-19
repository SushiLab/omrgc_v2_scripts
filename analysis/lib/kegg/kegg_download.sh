#!/bin/bash
echo "Downloading pathway list"
curl -# http://rest.kegg.jp/list/pathway > ./pathway.list

echo "Downloading module list"
curl -# http://rest.kegg.jp/list/module > ./module.list

echo "Downloading reaction list"
curl -# http://rest.kegg.jp/list/reaction > ./reaction.list

echo "Downloading ko list"
curl -# http://rest.kegg.jp/list/ko > ./ko.list

echo "Downloading ko-pathway link"
curl -# http://rest.kegg.jp/link/ko/pathway > ./ko_pathway.link

echo "Downloading ko-module link"
curl -# http://rest.kegg.jp/link/ko/module > ./ko_module.link

echo "Downloading ko-reaction link"
curl -# http://rest.kegg.jp/link/ko/reaction > ./ko_reaction.link

echo "Downloading module-pathway link"
curl -# http://rest.kegg.jp/link/module/pathway > ./module_pathway.link
