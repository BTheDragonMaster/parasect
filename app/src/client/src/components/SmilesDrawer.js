import React, { useEffect } from "react";
import { Box } from "@mui/material";
import SmilesDrawer from "smiles-drawer";

/**
 * component to draw a molecule from a SMILES string.
 * 
 * @param {number} identifier - Unique identifier for the component.
 * @param {string} smilesStr - SMILES string of the molecule.
 * @returns {React.ReactElement} - The component showing the molecule.
 */ 
const SmileDrawerContainer = ({ identifier, smilesStr, width, height }) => {
    // create a new drawer instance
    let drawer = new SmilesDrawer.SvgDrawer({ width: width, height: height });

    // draw the molecule when the component is mounted
    useEffect(() => {
        SmilesDrawer.parse(smilesStr, function (tree) {
            drawer.draw(tree, `structure-svg-${identifier}`, "light");
        });
    }, [smilesStr]); // re-draw the molecule when the SMILES string changes

    return (
        <Box key={identifier} sx={{ width: width, height: height }}>
            <svg id={`structure-svg-${identifier}`}/>
        </Box>
    );
};

export default SmileDrawerContainer;
