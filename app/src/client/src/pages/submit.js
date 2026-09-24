import React, { useRef, useState } from 'react';
import { toast } from 'react-toastify';
import { Box, Button, Divider, TextField, FormControl, FormLabel, RadioGroup, FormControlLabel, Radio, CircularProgress, Typography, Input } from '@mui/material';
import { MdSettings, MdBugReport } from 'react-icons/md';

import ExampleInputPicker from '../components/ExampleInputPicker';
import SettingsModal from '../components/SettingsModal';

const gbkExtensions = ['.gbk', '.gb', '.gbff', '.genbank'];
const fastaExtensions = ['.fasta', '.fa', '.faa'];

/**
 * Detect whether input is FASTA or GBK, from its content or else its file name.
 *
 * The backend parses input strictly as the selected type, so a GBK parsed as
 * FASTA yields no sequences and fails deep inside HMMER.
 *
 * @param {string} content - The input text.
 * @param {string} [fileName] - The name of the uploaded file, if any.
 * @returns {string|null} - 'fasta', 'gbk', or null if undetermined.
 */
const detectInputType = (content, fileName = '') => {
    const start = content.trimStart();
    if (start.startsWith('LOCUS')) return 'gbk';
    if (start.startsWith('>')) return 'fasta';

    const name = fileName.toLowerCase();
    if (gbkExtensions.some((ext) => name.endsWith(ext))) return 'gbk';
    if (fastaExtensions.some((ext) => name.endsWith(ext))) return 'fasta';

    return null;
};

/**
 * Component to submit data to the server.
 *
 * @param {Object} props - The props of the component.
 * @param {string} props.imageSrc - The path to the image.
 * @param {string} props.label - The label for the radio button.
 * @returns {React.ReactElement} - The submit component.
 */
const RadioLabel = ({ imageSrc, label }) => (
    <Box display="flex" alignItems="center">
        <Box
            component="img"
            src={imageSrc}
            alt=""
            sx={{ width: 40, height: 40, marginRight: 1, borderRadius: '6px' }}
        />
        <Typography variant="body1">{label}</Typography>
    </Box>
);


/**
 * Component to submit data to the server.
 * 
 * @returns {React.ReactElement} - The submit component.
 */
const Submit = () => {
    // page state
    const [isLoading, setIsLoading] = useState(false);

    // input method and type
    const [inputMethod, setInputMethod] = useState('upload'); // 'paste' or 'upload'
    const [selectedInputType, setSelectedInputType] = useState('fasta'); // fasta or gbk
    const [selectedInput, setSelectedInput] = useState('');
    const [loadedExampleFileName, setLoadedExampleFileName] = useState(null);
    const fileInputRef = useRef(null);

    // options
    const [selectedModel, setSelectedModel] = useState('parasAllSubstrates'); // parasAllSubstrates, parasCommonSubstrates, or parasect
    const [useStructureGuidedAlignment, setUseStructureGuidedAlignment] = useState(false);

    // SMILES file state (for PARASECT model)
    const [smilesFileContent, setSmilesFileContent] = useState(''); // Stores content of uploaded SMILES file
    const [useOnlyUploadedSubstrates, setUseOnlyUploadedSubstrates] = useState(false); // Checkbox state
    const [uploadedSubstratesFileContentHasHeader, setUploadedSubstratesFileContentHasHeader] = useState(true); // Checkbox state

    // modal state
    const [openSettingsModal, setOpenSettingsModal] = useState(false);

    // open and close modal handlers
    const handleOpenSettingsModal = () => setOpenSettingsModal(true);
    const handleCloseSettingsModal = () => setOpenSettingsModal(false);

    // load an example picked from the menu
    function handleLoadExample(example) {
        setSelectedInputType(example.inputType);
        setSelectedInput(example.content);
        setLoadedExampleFileName(example.fileName);
        if (fileInputRef.current) fileInputRef.current.value = '';
    };

    // refresh the page
    function handleRefresh() {
        // remove results from local storage
        localStorage.removeItem('results');

        // reload the page
        window.location.reload();
    };

    // handle file upload
    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (file) {
            setLoadedExampleFileName(null);
            const reader = new FileReader();
            reader.onload = function (event) {
                const fileContent = event.target.result;
                setSelectedInput(fileContent); // set the file content into selectedInput

                // match the input type to the file, so a GBK isn't parsed as FASTA
                const detectedType = detectInputType(fileContent, file.name);
                if (detectedType) setSelectedInputType(detectedType);
            };
            reader.readAsText(file);
        };
    };

    // handle form submission
    const handleSubmit = async () => {
        setIsLoading(true);

        const data = {
            selectedInputType: selectedInputType,
            selectedInput: selectedInput,
            selectedModel: selectedModel,
            useStructureGuidedAlignment: useStructureGuidedAlignment,
            smilesFileContent: smilesFileContent,
            useOnlyUploadedSubstrates: useOnlyUploadedSubstrates,
            uploadedSubstratesFileContentHasHeader: uploadedSubstratesFileContentHasHeader,
        };

        try {
            const response = await fetch('/api/submit_raw', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error('Network response was not ok!');
            };

            const json = await response.json();

            if (json.status === 'success') {
                const jobId = json.payload['jobId'];
                window.location.href = `/results/${jobId}`;
            } else if (json.status === 'warning') {
                toast.warn(json.message);
            } else if (json.status === 'failure') {
                toast.error(json.message);
            };
        } catch (error) {
            console.error('Error:', error);
            toast.error(error.message);
        };

        setIsLoading(false);
    };

    return (
        <>
            <Box
                display='flex'
                flexDirection='column'
                alignItems='left'
                sx={{ px: { xs: 2, sm: 4 }, py: 4, maxWidth: 840, width: '100%' }}
                margin='auto'
            >
                <Typography variant='h4' gutterBottom>
                    Submit your data
                </Typography>
                <Divider />
                
                {/* input method selection */}
                <FormControl component='fieldset' sx={{ mt: 3 }}>
                    <FormLabel component='legend'>Input method</FormLabel>
                    <RadioGroup
                        row
                        value={inputMethod}
                        onChange={(e) => setInputMethod(e.target.value)}
                    >
                        <FormControlLabel value='upload' control={<Radio />} label='Upload file' />
                        <FormControlLabel value='paste' control={<Radio />} label='Paste data' />
                    </RadioGroup>
                </FormControl>

                {/* input type selection */}
                <FormControl component='fieldset' margin='normal'>
                    <FormLabel component='legend'>Input type</FormLabel>
                    <RadioGroup
                        row
                        value={selectedInputType}
                        onChange={(e) => setSelectedInputType(e.target.value)}
                    >
                        <FormControlLabel value='fasta' control={<Radio />} label='FASTA' />
                        <FormControlLabel value='gbk' control={<Radio />} label='GBK' />
                    </RadioGroup>
                </FormControl>
                
                {/* text field or file upload */}
                <Box margin={1} >
                    {inputMethod === 'paste' ? (
                        <TextField
                            label='Input adenylation domain data (FASTA or GBK)'
                            multiline
                            rows={8}
                            fullWidth
                            variant='outlined'
                            value={selectedInput}
                            onChange={(e) => {
                                setSelectedInput(e.target.value);
                                setLoadedExampleFileName(null);
                                const detectedType = detectInputType(e.target.value);
                                if (detectedType) setSelectedInputType(detectedType);
                            }}
                            margin='normal'
                            placeholder='Paste your sequence here'
                        />
                    ) : (
                        <Box width='100%' sx={{ mt: 3, mb: 3 }}>
                            <Typography variant='body1' gutterBottom>
                                Upload your {selectedInputType.toUpperCase()} file:
                            </Typography>
                            <Input
                                type='file'
                                inputProps={{ accept: [...fastaExtensions, ...gbkExtensions].join(',') }} // accept FASTA or GBK files
                                onChange={handleFileUpload}
                                inputRef={fileInputRef}
                            />
                            {loadedExampleFileName && (
                                <Typography variant='body2' color='textSecondary' sx={{ mt: 1 }}>
                                    Using example file: <b>{loadedExampleFileName}</b>
                                </Typography>
                            )}
                        </Box>
                    )}

                    {/* load example button */}
                    <ExampleInputPicker onLoad={handleLoadExample} />
                </Box>

                {/* model selection */}
                <FormControl component='fieldset' margin='normal'>
                    <FormLabel component='legend'>Select model</FormLabel>
                    <RadioGroup
                        value={selectedModel}
                        onChange={(e) => setSelectedModel(e.target.value)}
                    >
                        <FormControlLabel 
                            value='parasAllSubstrates' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/paras.png'} 
                                    label={'PARAS (all substrates): predict adenylation domain substrate specificty for a pre-defined list of 223 substrates'}
                                />} 
                        />
                        <FormControlLabel 
                            value='parasCommonSubstrates' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/paras.png'} 
                                    label={'PARAS (common substrates): predict adenylation domain substrate specificty for a pre-defined list of 34 common substrates'}
                                />}
                        />
                        <FormControlLabel 
                            value='parasect' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/parasect.png'} 
                                    label={'PARASECT: predict if a list of pre-defined and/or user-supplied substrates interact with the adenylation domains. This model is trained on fungal and bacterial data.'}
                                />}
                        />
                        <FormControlLabel 
                            value='parasectBacterial' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/parasect.png'} 
                                    label={'PARASECT (bacterial): predict if a list of pre-defined and/or user-supplied substrates interact with the adenylation domains. This model is trained on bacterial data only.'}
                                />
                            }
                        />
                    </RadioGroup>
                </FormControl>

                {/* settings, submit, and refresh buttons */}
                <Box mt={4} display='flex' flexWrap='wrap' justifyContent='left' width='100%' gap={2}>
                    <Button
                        variant='contained'
                        color='primary'
                        startIcon={<MdSettings size={20} style={{ fill: 'currentColor' }}/>}
                        onClick={handleOpenSettingsModal}
                    >
                        Settings
                    </Button>
                    <Button 
                        variant='contained' 
                        color='primary' 
                        onClick={handleRefresh}
                    >
                        Refresh
                    </Button>
                    <Button
                        variant='contained'
                        color='secondary'
                        onClick={handleSubmit}
                        disabled={isLoading || !selectedInput}
                    >
                        {isLoading ? <CircularProgress size={24} /> : 'Submit'}
                    </Button>
                </Box>
            </Box>

            {/* settings modal */}
            <SettingsModal
                openSettingsModal={openSettingsModal}
                handleCloseSettingsModal={handleCloseSettingsModal}
                selectedModel={selectedModel}
                useStructureGuidedAlignment={useStructureGuidedAlignment}
                setUseStructureGuidedAlignment={setUseStructureGuidedAlignment}
                smilesFileContent={smilesFileContent}
                setSmilesFileContent={setSmilesFileContent}
                useOnlyUploadedSubstrates={useOnlyUploadedSubstrates}
                setUseOnlyUploadedSubstrates={setUseOnlyUploadedSubstrates}
                uploadedSubstratesFileContentHasHeader={uploadedSubstratesFileContentHasHeader}
                setUploadedSubstratesFileContentHasHeader={setUploadedSubstratesFileContentHasHeader}
            />
        </>
    );
};

export default Submit;
