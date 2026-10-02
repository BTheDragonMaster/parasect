import React, { useState } from 'react';
import { toast } from 'react-toastify';
import { Button, CircularProgress, ListItemText, Menu, MenuItem } from '@mui/material';
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';

/**
 * Button with a menu of example inputs served by the backend.
 *
 * @param {Object} props - The props of the component.
 * @param {Function} props.onLoad - Called with { id, label, fileName, inputType, content } of the picked example.
 * @returns {React.ReactElement} - The example picker.
 */
const ExampleInputPicker = ({ onLoad }) => {
    const [anchorEl, setAnchorEl] = useState(null);
    const [examples, setExamples] = useState(null);
    const [isLoading, setIsLoading] = useState(false);

    const fetchPayload = async (url) => {
        const response = await fetch(url);
        if (!response.ok) throw new Error('Network response was not ok!');

        const json = await response.json();
        if (json.status !== 'success') throw new Error(json.message);
        return json.payload;
    };

    const handleOpen = async (e) => {
        setAnchorEl(e.currentTarget);
        if (examples) return;

        setIsLoading(true);
        try {
            setExamples((await fetchPayload('/api/example_inputs')).examples);
        } catch (error) {
            toast.error(`Could not load the list of examples: ${error.message}`);
            setAnchorEl(null);
        };
        setIsLoading(false);
    };

    const handlePick = async (id) => {
        setAnchorEl(null);
        setIsLoading(true);
        try {
            onLoad(await fetchPayload(`/api/example_inputs/${encodeURIComponent(id)}`));
        } catch (error) {
            toast.error(`Could not load the example: ${error.message}`);
        };
        setIsLoading(false);
    };

    return (
        <>
            <Button
                variant='text'
                color='primary'
                onClick={handleOpen}
                disabled={isLoading}
                endIcon={isLoading ? <CircularProgress size={16} /> : <ArrowDropDownIcon />}
            >
                Load example input
            </Button>
            <Menu
                anchorEl={anchorEl}
                open={Boolean(anchorEl) && Boolean(examples)}
                onClose={() => setAnchorEl(null)}
            >
                {(examples || []).map(({ id, label, description }) => (
                    <MenuItem key={id} onClick={() => handlePick(id)} sx={{ maxWidth: 420, whiteSpace: 'normal' }}>
                        <ListItemText primary={label} secondary={description} />
                    </MenuItem>
                ))}
            </Menu>
        </>
    );
};

export default ExampleInputPicker;
