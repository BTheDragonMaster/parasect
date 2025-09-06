import React, { useEffect, useState, useMemo } from 'react';
import { useParams } from 'react-router-dom';
import { toast } from 'react-toastify';
import { Box, IconButton, Divider, Typography, Button, Modal, Tooltip, TextField, Stack, Chip, CircularProgress } from '@mui/material';
import { CheckCircle, ErrorOutline } from "@mui/icons-material";
import { MdClose } from 'react-icons/md';

import { useNavigate } from "react-router-dom";

import Loading from '../components/Loading';
import ProteinTile from '../components/ProteinTile';

import Turnstile from 'react-turnstile';


const SITE_KEY = process.env.REACT_APP_TURNSTILE_SITE_KEY;


function SubmitAnnotationsModal({ open, onClose, proteinAnnotations }) {
    const navigate = useNavigate();
    const [captchaToken, setCaptchaToken] = useState(null);
    const [submitting, setSubmitting] = useState(false);
    const [turnstileKey, setTurnstileKey] = useState(0); // force reset when modal re-opens

    // ORCID state
    const [orcidInput, setOrcidInput] = useState("");
    const [orcidNormalized, setOrcidNormalized] = useState("");
    const [orcidError, setOrcidError] = useState("");

    // Reference (DOI / PMID)
    const [refInput, setRefInput] = useState("");
    const [references, setReferences] = useState([]); // [{key,type:'doi'|'pmid', id, status:'pending'|'valid'|'invalid'|'error', title?, url?}]
    const MAX_REFS = 20;

    useEffect(() => {
        if (open) {
            // reset widget each time the modal opens
            setCaptchaToken(null);
            setTurnstileKey(k => k + 1);

            setOrcidInput("");
            setOrcidNormalized("");
            setOrcidError("");

            setRefInput("");
            setReferences([]);
        }
    }, [open]);

    // ----- ORCID helpers -----
    const normalizeOrcid = (value) => {
        if (!value) return '';
        // Strip URL prefix if present
        const cleaned = value
            .trim()
            .replace(/^https?:\/\/(www\.)?orcid\.org\//i, '')
            .replace(/[^0-9xX]/g, ''); // keep digits and X
        // Allow last char X, others must be digits; length should be 16
        return cleaned.toUpperCase();
    };

    const formatOrcidPretty = (digits16) => {
        // returns 0000-0000-0000-0000
        return `${digits16.slice(0,4)}-${digits16.slice(4,8)}-${digits16.slice(8,12)}-${digits16.slice(12,16)}`;
    };

    // ISO 7064 mod 11-2 checksum for ORCID (last char may be 'X' which means 10)
    const isValidOrcid = (digits16) => {
        if (!/^\d{15}[\dX]$/.test(digits16)) return false;
        let total = 0;
        for (let i = 0; i < 15; i++) {
            total = (total + parseInt(digits16[i], 10)) * 2;
        }
        const remainder = total % 11;
        const result = (12 - remainder) % 11;
        const checkChar = result === 10 ? 'X' : String(result);
        return checkChar === digits16[15];
    };

    const handleOrcidChange = (e) => {
        const raw = e.target.value;
        setOrcidInput(raw);

        const normalized = normalizeOrcid(raw);
        if (normalized.length === 16 && isValidOrcid(normalized)) {
            setOrcidNormalized(formatOrcidPretty(normalized));
            setOrcidError('');
        } else {
            setOrcidNormalized('');
            // only show error once user has typed enough chars that could plausibly be an ORCID
            if (normalized.length >= 12) {
                setOrcidError('Invalid ORCID. Expected format 0000-0000-0000-000X with a valid checksum.');
            } else {
                setOrcidError('');
            }
        }
    };

    // ----- Reference helpers -----
    const parseRefToken = (token) => {
        const t = token.trim();
        if (!t) return null;

        // Strip URL wrappers
        const stripped = t
            .replace(/^https?:\/\/(dx\.)?doi\.org\//i, "")
            .replace(/^https?:\/\/pubmed\.ncbi\.nlm\.nih\.gov\//i, "")
            .replace(/^pmid:\s*/i, "");

        // DOI regex (from Crossref guide)
        const doiRegex = /^10\.\d{4,9}\/[-._;()/:A-Z0-9]+$/i;
        if (doiRegex.test(stripped)) {
            return { type: "doi", id: stripped };
        }

        // PMID (digits only, typical up to 8-9 digits)
        const pmidRegex = /^\d{1,9}$/;
        if (pmidRegex.test(stripped)) {
            return { type: "pmid", id: stripped };
        }

        return null;
    };

    const addReferencesFromInput = async () => {
        if (!refInput.trim()) return;

        const tokens = refInput
            .split(/[\s,;\n]+/) // spaces, commas, semicolons, newlines
            .map((x) => x.trim())
            .filter(Boolean);

        let current = [...references];
        const seen = new Set(current.map((r) => `${r.type}:${r.id}`));

        for (const tok of tokens) {
            if (current.length >= MAX_REFS) break;
            const parsed = parseRefToken(tok);
            if (!parsed) continue;
            const key = `${parsed.type}:${parsed.id}`.toLowerCase();
            if (seen.has(key)) continue;

            seen.add(key);
            current.push({
                key,
                ...parsed,
                status: "pending",
                title: undefined,
                url: undefined,
            });
        }

        if (current.length > references.length) {
            setReferences(current);
            setRefInput("");
            // validate the new ones
            const newOnes = current.filter((r) => r.status === "pending");
            for (const ref of newOnes) {
                validateReference(ref);
            }
        } else {
            // nothing added
            setRefInput("");
        }
    };

    const removeReference = (key) => {
        setReferences((prev) => prev.filter((r) => r.key !== key));
    };

    const markRef = (key, patch) => {
        setReferences((prev) => prev.map((r) => (r.key === key ? { ...r, ...patch } : r)));
    };

    // Crossref for DOI
    const validateViaCrossref = async (doi) => {
        const url = `https://api.crossref.org/works/${encodeURIComponent(doi)}`;
        const res = await fetch(url, { headers: { Accept: "application/json" } });
        if (!res.ok) throw new Error(`Crossref ${res.status}`);
        const json = await res.json();
        const msg = json?.message;
        if (!msg) return { valid: false };

        const type = msg.type; // 'journal-article', 'posted-content' (preprint), etc.
        const issued = msg.issued; // date parts if available
        const isPreprint = type === "posted-content";
        const isPublished = !!issued && !isPreprint;

        const title = (Array.isArray(msg.title) && msg.title[0]) || doi;
        const urlOut = msg.URL || `https://doi.org/${doi}`;

        return { valid: Boolean(isPublished), title, url: urlOut };
    };

    // PubMed E-utilities for PMID
    const validateViaPubMed = async (pmid) => {
        const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=${encodeURIComponent(
            pmid
        )}&retmode=json`;
        const res = await fetch(url);
        if (!res.ok) throw new Error(`NCBI ${res.status}`);
        const json = await res.json();
        const rec = json?.result?.[pmid];
        if (!rec) return { valid: false };

        const title = rec.title || `PMID:${pmid}`;
        // consider "published" if has a pubdate and is not preprint. PubMed labels preprints differently
        const hasDate = Boolean(rec.pubdate || rec.epubdate || rec.sortpubdate);
        const isPreprint =
            Array.isArray(rec.pubtype) && rec.pubtype.some((pt) => /preprint/i.test(pt)); // conservative
        const valid = hasDate && !isPreprint;

        const urlOut = `https://pubmed.ncbi.nlm.nih.gov/${pmid}/`;
        return { valid, title, url: urlOut };
    };

    const validateReference = async (ref) => {
        try {
            if (ref.type === "doi") {
                const { valid, title, url } = await validateViaCrossref(ref.id);
                markRef(ref.key, { status: valid ? "valid" : "invalid", title, url });
            } else if (ref.type === "pmid") {
                const { valid, title, url } = await validateViaPubMed(ref.id);
                markRef(ref.key, { status: valid ? "valid" : "invalid", title, url });
            } else {
                markRef(ref.key, { status: "invalid" });
            }
        } catch (e) {
            // Fallback to backend if CORS/network fails
            // TODO: implement server endpint to validate reference
            try {
                const res = await fetch("/api/validate_reference", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({ type: ref.type, id: ref.id }),
                });
                if (res.ok) {
                    const j = await res.json();
                    const { valid, title, url } = j || {};
                    markRef(ref.key, { status: valid ? "valid" : "invalid", title, url });
                } else {
                    markRef(ref.key, { status: "error" });
                }
            } catch {
                markRef(ref.key, { status: "error" });
            }
        }
    };
    
    const refsPending = references.some((r) => r.status === "pending");
    const refsInvalid = references.some((r) => r.status === "invalid" || r.status === "error");
    const refsValidPayload = useMemo(
        () =>
            references
                .filter((r) => r.status === "valid")
                .map(({ type, id, title, url }) => ({ type, id, title, url })),
        [references]
    );

    // ORCID is optional, need at least one valid reference. However if an ORCID is provided, it must be valid.
    const canSubmit =
        !!captchaToken &&
        !submitting &&
        Object.keys(proteinAnnotations).length > 0 &&
        (!orcidInput.trim() || (!!orcidNormalized && !orcidError)) &&
        !refsPending &&
        !refsInvalid &&
        refsValidPayload.length > 0;

    // Submit updated protein data 
    const handleSubmit = async () => {
        if (submitting) return; // prevent multiple submissions

        if (Object.keys(proteinAnnotations).length === 0) {
            toast.error("No annotations to submit");
            return;
        };

        if (!captchaToken) return;

        // Require valid ORCID
        if (!orcidNormalized || orcidError) {
            setOrcidError('Please provide a valid ORCID.');
            return;
        }

        // Require valid references that are not pending
        if (refsPending || refsInvalid) {
            toast.error("Please resolve or remove invalid/pending references.");
            return;
        }

        setSubmitting(true);

        let prUrl = null;
    
        try {
            const res = await fetch("/api/submit_annotations", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({
                    annotations: proteinAnnotations,
                    turnstileToken: captchaToken,
                    orcid: orcidNormalized,
                    references: refsValidPayload
                }),
                credentials: 'include' // include cookies
            });

            const data = await res.json();

            if (!res.ok) {
                throw new Error(data?.error || "Submission failed")
            }

            prUrl = data?.pr_url || null;
            console.log(prUrl);

            // open PR url in separate tab
            if (prUrl) {
                window.open(prUrl, "_blank");
            }

            // success UI
            toast.success("Annotations submitted successfully");
            onClose?.();
            navigate("/data_annotation");
        } catch(e) {
            console.error(e);
            toast.error(`Error: ${e.message}`, { autoClose: false });
        } finally {
            setSubmitting(false);
            setCaptchaToken(null);
            setTurnstileKey(k => k + 1);  // reset widget after submit
        }
    };

    return (
        <Modal open={open} onClose={onClose}>
            <Box
                width={800}
                bgcolor='white.main'
                mx='auto'
                my={10}
                borderRadius={4}
                boxShadow={3}
            >
                <Box
                    sx={{
                        backgroundColor: 'secondary.main',
                        borderTopLeftRadius: '14px',
                        borderTopRightRadius: '14px',
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        padding: 1,
                    }}
                >
                    <Typography 
                        variant='h5' 
                        gutterBottom
                        sx={{ 
                            color: 'black.main', 
                            textAlign: 'center',
                            pl: 2,
                            pt: 2, 
                        }}
                    >
                        Submit annotations
                    </Typography>
                    <IconButton onClick={onClose}>
                        <MdClose size={24} />
                    </IconButton>
                </Box>

                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        alignItems: 'flex-start',
                        gap: 3,
                        p: 4,
                    }}
                >
                    {/* ORCID field */}
                    <Box sx={{ width: "100%" }}> 
                        <Typography variant="subtitle1" sx={{ mb: 1}}>
                            Contributor ORCID - optional
                        </Typography>
                        <TextField 
                            label="Contributor ORCID"
                            placeholder="e.g., 0000-0002-1825-0097 or https://orcid.org/0000-0002-1825-0097"
                            fullWidth
                            value={orcidInput}
                            onChange={handleOrcidChange}
                            error={!!orcidError}
                            helperText={orcidError || (orcidNormalized ? `Using: ${orcidNormalized}` : 'Enter your 16-character ORCID')}
                        />
                    </Box>

                    {/* References */}
                    <Box sx={{ width: "100%" }}> 
                        <Typography variant="subtitle1" sx={{ mb: 1}}>
                            Related publications (DOI or PubMed ID) - at least one required
                        </Typography>

                        <Box sx={{ display: "flex", gap: 1 }}>
                            <TextField 
                                fullWidth
                                label={`Enter DOI or PMID${references.length ? " (you can paste multiple)" : ""}`}
                                placeholder="10.1038/nature14539, 12345678, https://doi.org/10.1093/nar/gkv123, PMID: 9876543"
                                value={refInput}
                                onChange={(e) => setRefInput(e.target.value)}
                                helperText={`${references.length}/${MAX_REFS} added`}
                            />
                            <Button
                                variant="outlined"
                                onClick={addReferencesFromInput}
                                disabled={!refInput.trim() || references.length >= MAX_REFS}
                                sx={{ whiteSpace: "nowrap", height: '56px' }}
                            >
                                Add
                            </Button>
                        </Box>

                        {!!references.length && (
                            <Stack direction="row" spacing={1} useFlexGap flexWrap="wrap" sx={{ mt: 1 }}>
                                {references.map((r) => {
                                    const color =
                                        r.status === "valid" ? "success" : r.status === "pending" ? "default" : "error";
                                    const icon =
                                        r.status === "valid" ? (
                                            <CheckCircle fontSize="small" />
                                        ) : r.status === "pending" ? (
                                            <CircularProgress size={14} />
                                        ) : (
                                            <ErrorOutline fontSize="small" />
                                        );
                                    const label =
                                        r.type === "doi" ? `DOI:${r.id}` : r.type === "pmid" ? `PMID:${r.id}` : r.id;

                                    return (
                                        <Tooltip
                                            key={r.key}
                                            title={r.title ? `${r.title}${r.url ? `\n${r.url}` : ""}` : ""}
                                            arrow
                                        >
                                        <Chip
                                            icon={icon}
                                            label={label}
                                            color={color}
                                            onDelete={() => removeReference(r.key)}
                                            variant={r.status === "valid" ? "filled" : "outlined"}
                                        />
                                        </Tooltip>
                                    );
                                })}
                            </Stack>
                        )}

                        {refsInvalid && (
                            <Typography variant="body2" sx={{ mt: 1 }} color="error">
                                Some references are invalid or not recognized as published articles. Remove or correct them to submit.
                            </Typography>
                        )}
                        {refsPending && (
                            <Typography variant="body2" sx={{ mt: 1 }}>
                                Validating references…
                            </Typography>
                        )}
                    </Box>

                    {/* Turnstile widget */}
                    <Turnstile
                        key={turnstileKey}
                        sitekey={SITE_KEY}
                        onVerify={(token) => setCaptchaToken(token)}
                        onExpire={() => setCaptchaToken(null)}
                        onError={() => setCaptchaToken(null)}
                        options={{ theme: 'auto', appearance: 'always' }}
                    />
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={handleSubmit}
                        disabled={!canSubmit}
                    >
                        {submitting 
                            ? 'Submitting…' 
                            : `Submit ${Object.keys(proteinAnnotations).length} annotated protein(s)`}
                    </Button>
                </Box>
            </Box>
        </Modal>
    )  
};


/**
 * Component to display the results of the prediction.
 *
 * @returns {React.ReactElement} - The results component.
 */

const AnnotationEditor = () => {
    // get job ID from URL
    const { jobId } = useParams();

    // state to store results
    const [results, setResults] = useState(null);

    // state to keep track of loading state
    const [isLoading, setIsLoading] = useState(true);

    const [proteinAnnotations, setProteinAnnotations] = useState({});

    // submission modal state
    const [openAnnotationsSubmissionModal, setOpenAnnotationsSubmissionModal] = useState(false);

    const handleCloseAnnotationsSubmissionModal = () => {
        setOpenAnnotationsSubmissionModal(false);
    };

    {/* For collecting protein annotations */}

    const handleProteinAnnotationChange = (proteinId, data) => {
    setProteinAnnotations((prev) => {
        const updated = { ...prev };

        // Check if domains object exists and has keys
        const hasDomainAnnotations = data.domains && Object.keys(data.domains).length > 0;

        if (data && hasDomainAnnotations) {
            updated[proteinId] = data;
        } else {
            // Remove if no domain annotations (ignore synonym)
            delete updated[proteinId];
        }

        return updated;
    });
};

    const OpenAnnotationsSubmissionsModal = () => {
        setOpenAnnotationsSubmissionModal(true);
    };

    // fetch results from local storage
    useEffect(() => {
        let intervalId;

        const fetchResult = async () => {
            try {
                const response = await fetch(`/api/retrieve/${jobId}`);
                if (!response.ok) {
                    throw new Error('failed to fetch results');
                }

                const data = await response.json();

                if (data.status === 'success') {
                    const results = data['payload']['results']

                    // set states
                    setResults(results);
                    setIsLoading(false);
                    clearInterval(intervalId);
                } else if (data.status === 'failure') {
                    throw new Error(data.message);
                } // else keep polling
            } catch (error) {
                toast.error(
                    <>
                      {error.message}<br /><br />
                      If you feel this is an error, or if you need assistance, please contact the developers in GitHub issues by selecting 'Report an issue' in the app bar at the top left of this page and posting your issue or question.
                    </>,
                    { autoClose: false }
                );
                setIsLoading(false);
                clearInterval(intervalId);
            }
        };

        if (jobId) {
            // poll every second (1000 milliseconds)
            intervalId = setInterval(fetchResult, 1000);
        }

        // clear interval when component unmounts
        return () => clearInterval(intervalId);

    }, [jobId]);

    // render loading spinner while fetching results
    if (isLoading) {
        return (
            <Box
                display='flex'
                flexDirection='column'
                justifyContent='center'
                alignItems='center'
                minHeight='80vh'
            >
                <Loading
                    frame1='paras_loading_1.png'
                    frame2='paras_loading_2.png'
                />
                <p>Extracting A-domains and making PARAS predictions...</p>
            </Box>
        );
    }

    // render message if no results are found
    if (!results) {
        return (
            <Box
                display='flex'
                justifyContent='center'
                alignItems='center'
                minHeight='80vh'
            >
                <p>No results found for job ID {jobId}</p>
            </Box>
        );
    }

    // render results if available
    return (
        <>
            <Box
                display='flex'
                flexDirection='column'
                overflow='hidden'
            >
                <Box
                    display='flex'
                    flexDirection='column'
                    alignItems='left'
                    margin={4}
                >
                    <Box sx={{ display: 'flex', flexDirection: 'row', gap: 1 }}>
                        <Typography variant='h4' gutterBottom>
                            Extracted adenylation domains
                        </Typography>
                    </Box>
                    <Divider />

                    <Box sx={{ mt: 4 }}>
                        <Typography variant='body1' gutterBottom>
                            In total, {results.length} of the submitted proteins contain adenylation domains. The tiles below show the domains for each protein. You can scroll horizontally to view all proteins.
                        </Typography>

                        <Typography variant='body1' gutterBottom>
                            Domains which already exist in the PARAS/PARASECT dataset are displayed in grey. New domains are displayed in yellow. Please review the new domains and provide annotations where possible. You can proceed with submitting domains once they are annotated.
                        </Typography>

                        <Typography variant='body1' gutterBottom>
                        </Typography>
                        <Box sx={{ mt: 4 }}>
                            <Tooltip title={Object.keys(proteinAnnotations).length === 0 ? "Nothing to submit! Have you annotated all new domains?" : ""} arrow>
                                {/* Wrap Button in a span to support tooltip for disabled state */}
                                <span>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={OpenAnnotationsSubmissionsModal}
                                        disabled={Object.keys(proteinAnnotations).length === 0}
                                    >
                                        Proceed with submission
                                    </Button>
                                </span>
                            </Tooltip>
                        </Box>
                    </Box>
                </Box>

                {/* display results in a row, one item per protein */}
                <Box
                    sx={{
                        overflowY: 'auto',
                        overflowX: 'auto',
                        backgroundColor: 'white.main',
                        flexDirection: 'row',
                        display: 'flex',
                        gap: '20px',
                        paddingLeft: '30px',
                        paddingRight: '30px',
                        paddingBottom: '20px',

                        // always show scrollbar
                        '&::-webkit-scrollbar': {
                            display: 'block',
                        },

                        // scrollbar style
                        '&::-webkit-scrollbar-thumb': {
                            backgroundColor: '#ceccca',
                            borderRadius: '10px',
                        },
                    }}
                >
                    {results.map((result, index) => (
                        <ProteinTile
                            key={index}
                            proteinResult={result}
                            onUpdateAnnotation={handleProteinAnnotationChange} />
                    ))}
                </Box>
            </Box>

            {/* Dialogue modal */}
            <SubmitAnnotationsModal 
                open={openAnnotationsSubmissionModal}
                onClose={handleCloseAnnotationsSubmissionModal}
                proteinAnnotations={proteinAnnotations}
            />
        </>
        
    );
};

export default AnnotationEditor;