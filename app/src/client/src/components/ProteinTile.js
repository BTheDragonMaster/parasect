import React, { useState } from 'react';
import { toast } from 'react-toastify';
import { Box, Button, MenuItem, Select, FormControl } from '@mui/material';
import { FaFingerprint, FaCopy } from 'react-icons/fa';

import SmileDrawerContainer from './SmilesDrawer';
import DomainTile from '../components/DomainTile';

/**
 * Component to display the results of the prediction.
 *
 * @param {Object} props - The component props.
 * @param {Object} props.result - The result object.
 * @returns {React.ReactElement} - The result tile component.
 */
const ProteinTile = ({ proteinResult }) => {

    return (
        <Box
            sx={{
                minWidth: '800px',
                maxWidth: '800px',
                borderRadius: '11px',
                boxShadow: '0px 4px 10px rgba(100, 84, 31, 0.5)',
            }}
        >

            {/* header with domain name and location */}
            <Box
                sx={{
                    backgroundColor: 'secondary.main',
                    color: 'black.main',
                    padding: '14px 8px',
                    display: 'flex',
                    marginBottom: 1,
                    borderTopLeftRadius: '10px',
                    borderTopRightRadius: '10px',
                    justifyContent: 'center',
                }}
            >
                {`${proteinResult['protein_name']}`}
            </Box>
            {/* Domains in the protein > 0 */}

            {proteinResult['results'].map((result, index) => (
                    <DomainTile key={index} result={result} />
                ))}


            }

        </Box>
    );
};

export default ProteinTile;
