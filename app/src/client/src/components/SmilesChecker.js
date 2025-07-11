import React, { useState } from 'react';
import {
  TextField,
  Box,
  CircularProgress,
  Typography,
  Alert,
  Button
} from '@mui/material';
import SmileDrawerContainer from './SmilesDrawer';

const SmilesChecker = ({ i, customSmiles, updateSubstrateField }) => {
  const [candidateCustomSmiles, setCandidateCustomSmiles] = useState(customSmiles || '');
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [invalidSmiles, setInvalidSmiles] = useState(false);
  const [lastCheckedSmiles, setLastCheckedSmiles] = useState('');

  const handleCheckSmiles = async () => {
  if (!candidateCustomSmiles) {
    setResult(null);
    setInvalidSmiles(false);
    setLastCheckedSmiles('');
    return;
  }

  setLoading(true);
  setError(null);

  try {
    const res = await fetch('/api/check_smiles', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: candidateCustomSmiles }),
    });

    if (!res.ok) throw new Error('Server error');

    const data = await res.json();

    if (data.invalid) {
      setInvalidSmiles(true);
      setResult(null);
      setLastCheckedSmiles('');
      return;
    }

    setInvalidSmiles(false);
    setResult(data);

    // Update the real SMILES value
    updateSubstrateField(i, 'customSmiles', candidateCustomSmiles);

    if (Array.isArray(data) && data.length > 0) {
      updateSubstrateField(i, 'substrateName', null);
      updateSubstrateField(i, 'substrateSmiles', null);
    }

    setLastCheckedSmiles(candidateCustomSmiles);  // ✅ Set the last checked SMILES
  } catch (err) {
    console.error(err);
    setResult(null);
    setError('An error occurred while checking the SMILES.');
    setLastCheckedSmiles('');
  } finally {
    setLoading(false);
  }
};

  const hasChanged = candidateCustomSmiles !== customSmiles;

  return (
    <Box sx={{ mt: 1 }}>
      <TextField
        label="Enter SMILES"
        fullWidth
        value={candidateCustomSmiles}
        onChange={(e) => setCandidateCustomSmiles(e.target.value)}
      />

      <Button
        variant="contained"
        sx={{ mt: 2 }}
        onClick={handleCheckSmiles}
        disabled={loading || !candidateCustomSmiles || !hasChanged}
      >
        {loading ? 'Checking...' : 'Check SMILES'}
      </Button>

      {loading && (
        <Box sx={{ mt: 1, display: 'flex', alignItems: 'center' }}>
          <CircularProgress size={20} />
          <Typography variant="caption" sx={{ ml: 1 }}>
            Checking...
          </Typography>
        </Box>
      )}

      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}

      {invalidSmiles && (
        <Alert severity="error" sx={{ mt: 2 }}>
          Invalid SMILES string.
        </Alert>
      )}

      {!loading && result && result.length > 0 && (
        <>
          <Alert severity="info" sx={{ mt: 2 }}>
            Molecule may already be in the dataset. Untick 'New substrate' and select the substrate from the dropdown menu.
          </Alert>

          <Box sx={{ mt: 2 }}>
            {result.map((substrate, idx) => (
              <Box
                key={idx}
                sx={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'space-between',
                  p: 2,
                  mb: 2,
                  border: '1px solid #ccc',
                  borderRadius: 2,
                  backgroundColor: '#f9f9f9',
                }}
              >
                <Box sx={{ flexGrow: 1 }}>
                  <Typography variant="subtitle2">
                    {substrate.name}
                  </Typography>
                </Box>
                <Box sx={{ ml: 2 }}>
                  <SmileDrawerContainer
                    identifier={`${i}-${idx}`}
                    smilesStr={substrate.smiles}
                    width={120}
                    height={100}
                  />
                </Box>
              </Box>
            ))}
          </Box>
        </>
      )}

      {!loading &&
            result &&
            result.length === 0 &&
            lastCheckedSmiles === candidateCustomSmiles &&
            !invalidSmiles && (
              <Typography variant="body2" sx={{ mt: 1 }}>
                Valid SMILES string
              </Typography>
          )}
    </Box>
  );
};

export default SmilesChecker;