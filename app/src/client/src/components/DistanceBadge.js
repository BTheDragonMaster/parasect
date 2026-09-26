import React from 'react';
import { Box, Tooltip, Typography } from '@mui/material';

import { useColorMode } from '../theme/ColorModeContext';
import { DISTANCE_RAMP } from '../theme';

/** Residues in an extended signature; a distance is a count out of this. */
export const SIGNATURE_LENGTH = 34;

/**
 * The five bands a Hamming distance falls into, quietest first.
 *
 * The cuts are about what the number means for a substrate call rather than
 * even fifths of 0-34: a handful of mismatches in 34 residues still leaves two
 * domains in the same specificity pocket, while past ten they are only
 * distantly related and nothing about the neighbour transfers. upTo is
 * inclusive, and the last band is open-ended.
 */
export const DISTANCE_BANDS = [
    { upTo: 0, label: 'identical', range: '0' },
    { upTo: 2, label: 'near-identical', range: '1-2' },
    { upTo: 5, label: 'close', range: '3-5' },
    { upTo: 10, label: 'related', range: '6-10' },
    { upTo: SIGNATURE_LENGTH, label: 'distant', range: '11+' },
];

/**
 * Band index for a distance, clamped to the ramp.
 *
 * @param {number} distance - mismatching residues, 0-34.
 * @returns {number} - index into DISTANCE_BANDS / DISTANCE_RAMP.
 */
export function distanceBand(distance) {
    const index = DISTANCE_BANDS.findIndex((band) => distance <= band.upTo);
    return index === -1 ? DISTANCE_BANDS.length - 1 : index;
}

/** Ramp step for a distance in the active colour mode. */
function useStep(distance) {
    const { mode } = useColorMode();
    return DISTANCE_RAMP[mode][distanceBand(distance)];
}

/**
 * A distance as a colour-ramped chip with the count inside it.
 *
 * The number is the point, the fill only says "how far" at a glance, and
 * carrying both means the badge still reads when the hue does not (colour
 * blindness, a printed figure, a greyscale screenshot).
 *
 * @param {object} props - component props.
 * @param {number} props.distance - mismatching residues, 0-34.
 * @returns {React.ReactElement} - the badge.
 */
export const DistanceBadge = ({ distance }) => {
    const step = useStep(distance);
    const band = DISTANCE_BANDS[distanceBand(distance)];
    const mismatches = `${distance} mismatch${distance === 1 ? '' : 'es'}`;

    return (
        <Tooltip title={`${mismatches} of ${SIGNATURE_LENGTH} residues - ${band.label}`}>
            <Box
                aria-label={`${mismatches}, ${band.label}`}
                sx={{
                    flexShrink: 0,
                    minWidth: 26,
                    height: 22,
                    px: 0.75,
                    borderRadius: '6px',
                    display: 'inline-flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    backgroundColor: step.fill,
                    color: step.text,
                    fontSize: '0.75rem',
                    fontWeight: 700,
                    fontVariantNumeric: 'tabular-nums',
                    lineHeight: 1,
                }}
            >
                {distance}
            </Box>
        </Tooltip>
    );
};

/**
 * Key for the distance ramp: the five bands as one strip, quietest first.
 *
 * Shown once above a list of badges rather than repeated per row. The badges
 * are the data, this only says what their fills mean.
 *
 * @returns {React.ReactElement} - the legend.
 */
export const DistanceLegend = () => {
    const { mode } = useColorMode();
    const ramp = DISTANCE_RAMP[mode];

    return (
        <Box>
            <Box sx={{ display: 'flex', gap: '2px' }}>
                {DISTANCE_BANDS.map((band, index) => (
                    <Tooltip key={band.label} title={`${band.range} of ${SIGNATURE_LENGTH} residues differ`}>
                        <Box
                            sx={{
                                flex: 1,
                                minWidth: 0,
                                py: 0.25,
                                borderRadius: index === 0 ? '4px 0 0 4px'
                                    : (index === DISTANCE_BANDS.length - 1 ? '0 4px 4px 0' : 0),
                                backgroundColor: ramp[index].fill,
                                color: ramp[index].text,
                                fontSize: '0.65rem',
                                fontWeight: 700,
                                fontVariantNumeric: 'tabular-nums',
                                textAlign: 'center',
                            }}
                        >
                            {band.range}
                        </Box>
                    </Tooltip>
                ))}
            </Box>
            <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
                {/* the ends are named; the three middle bands are read off the
                    ramp between them, which is what a ramp is for */}
                <Typography variant='caption' color='textSecondary'>
                    {DISTANCE_BANDS[0].label}
                </Typography>
                <Typography variant='caption' color='textSecondary'>
                    {DISTANCE_BANDS[DISTANCE_BANDS.length - 1].label}
                </Typography>
            </Box>
        </Box>
    );
};

export default DistanceBadge;
