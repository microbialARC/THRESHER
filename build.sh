#!/usr/bin/env bash

mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/workflow"
cp -r "${RECIPE_DIR}/thresher/workflow/" "${PREFIX}/workflow/"
cp "${RECIPE_DIR}/thresher/bin/cli.py" "${PREFIX}/bin/cli.py"
chmod +x "${PREFIX}/workflow"
chmod +x "${PREFIX}/bin/cli.py"