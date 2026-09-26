import React, { useEffect } from "react";
import { Box } from "@mui/material";
import SmilesDrawer from "smiles-drawer";

import { useColorMode } from "../theme/ColorModeContext";

/**
 * component to draw a molecule from a SMILES string.
 *
 * smiles-drawer paints carbons and bonds in the theme's ink, near-black under
 * its "light" theme, so the structure has to follow the app's colour mode or
 * it disappears into a dark card, leaving only the coloured heteroatoms.
 *
 * @param {number} identifier - Unique identifier for the component.
 * @param {string} smilesStr - SMILES string of the molecule.
 * @param {number} width - Drawing width in px.
 * @param {number} height - Drawing height in px.
 * @returns {React.ReactElement} - The component showing the molecule.
 */
const SmileDrawerContainer = ({ identifier, smilesStr, width, height }) => {
    const { mode } = useColorMode();

    // re-draw when the SMILES changes or the colour mode flips
    useEffect(() => {
        const target = document.getElementById(`structure-svg-${identifier}`);
        if (!target) return;
        // draw() appends into the element, so clear it or the old structure
        // stays underneath the new one on a redraw
        target.innerHTML = '';

        const drawer = new SmilesDrawer.SvgDrawer({ width, height });
        SmilesDrawer.parse(smilesStr, (tree) => {
            drawer.draw(tree, `structure-svg-${identifier}`, mode);
        });
    }, [smilesStr, identifier, width, height, mode]);

    return (
        <Box key={identifier} sx={{ width: width, height: height }}>
            <svg id={`structure-svg-${identifier}`}/>
        </Box>
    );
};

export default SmileDrawerContainer;
