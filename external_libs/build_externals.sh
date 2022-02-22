#!/usr/bin/env bash

./build_hypre.sh
./build_silo.sh

rm -r ./silo-4.10.2-bsd
rm -r ./hypre-2.11.2
