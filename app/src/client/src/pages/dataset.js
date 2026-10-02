import React, { useEffect, useMemo, useState } from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { toast } from 'react-toastify';
import {
    Box, Button, Card, CircularProgress, Divider, IconButton, Link, Stack, Tooltip, Typography,
} from '@mui/material';
import { useTheme } from '@mui/material/styles';
import DownloadIcon from '@mui/icons-material/Download';
import ContentCopyIcon from '@mui/icons-material/ContentCopy';
import StorageIcon from '@mui/icons-material/Storage';
import ExitIcon from '@mui/icons-material/ExitToApp';

import { DATABASE_LICENSE_URL, PUBLICATION_URL } from '../links';
import { CATEGORICAL } from '../theme';
import { schemaToSvg } from '../utils/erDiagram';
import { svgToPng } from '../utils/networkExport';
import { downloadBlob } from '../utils/zip';


const EXAMPLE_QUERY = `-- the ten substrates with the most annotated A domains
SELECT s.name, COUNT(*) AS domains
FROM substrate AS s
JOIN substrate_domain_association AS sda ON sda.substrate_name = s.name
GROUP BY s.name
ORDER BY domains DESC
LIMIT 10;`;

const formatBytes = (bytes) => (bytes >= 1024 * 1024
    ? `${(bytes / (1024 * 1024)).toFixed(1)} MB`
    : `${Math.ceil(bytes / 1024)} KB`);

/** Copy text to the clipboard and say so. */
const copyText = (text, what) => {
    navigator.clipboard.writeText(text)
        .then(() => toast.success(`Copied ${what} to clipboard`))
        .catch(() => toast.error(`Could not copy ${what}`));
};

/**
 * Page for downloading the full reference database, with a diagram of its schema.
 *
 * @returns {React.ReactElement} - The dataset download page.
 */
const Dataset = () => {
    const theme = useTheme();
    const mode = theme.palette.mode;
    const [info, setInfo] = useState(null);
    const [error, setError] = useState(null);

    useEffect(() => {
        fetch('/api/dataset/info')
            .then((response) => {
                if (!response.ok) throw new Error('Network response was not ok!');
                return response.json();
            })
            .then((json) => {
                if (json.status !== 'success') throw new Error(json.message);
                setInfo(json.payload);
            })
            .catch((err) => setError(err.message));
    }, []);

    // redrawn per colour mode; the download is whatever is on screen, as on the network page
    const diagram = useMemo(() => {
        if (!info) return null;
        return schemaToSvg({
            tables: info.tables,
            caption: {
                title: 'PARAS/PARASECT database schema',
                subtitle: `${info.fileName} / ${info.tables.length} tables`,
            },
            colors: {
                background: theme.palette.background.paper,
                paper: theme.palette.background.paper,
                header: theme.palette.surface.sunken,
                border: theme.palette.surface.borderStrong,
                text: theme.palette.text.primary,
                textSecondary: theme.palette.text.secondary,
                edge: theme.palette.text.secondary,
                primaryKey: theme.palette.primary.main,
                foreignKey: CATEGORICAL[mode][2],
            },
        });
    }, [info, theme, mode]);

    const baseName = info ? info.fileName.replace(/\.db$/, '') : 'parasect';

    const downloadSvg = () => {
        downloadBlob(new Blob([diagram.svg], { type: 'image/svg+xml' }), `${baseName}-schema.svg`);
    };

    const downloadPng = async () => {
        try {
            downloadBlob(await svgToPng(diagram.svg), `${baseName}-schema.png`);
        } catch (err) {
            toast.error(`Could not create the PNG: ${err.message}`);
        }
    };

    return (
        <Box sx={{ maxWidth: 1180, mx: 'auto', px: { xs: 2, sm: 4 }, py: 4 }}>
            <Typography variant='h4' gutterBottom>
                Download the dataset
            </Typography>
            <Divider />

            <Typography variant='body1' sx={{ mt: 3, maxWidth: 760 }}>
                The complete PARAS and PARASECT reference dataset is a single SQLite database: every adenylation
                domain with its signatures and substrate annotations, the proteins the domains come from, and their
                taxonomy. It is the same data the <Link component={RouterLink} to='/query_database'>query page</Link> searches. Open it with
                the <code>sqlite3</code> command-line tool, DB Browser for SQLite, or any SQLite library in Python or R.
            </Typography>

            {error && (
                <Typography color='error' sx={{ mt: 3 }}>
                    Could not load the dataset details: {error}
                </Typography>
            )}
            {!info && !error && <CircularProgress sx={{ mt: 4 }} />}

            {info && (
                <>
                    <Card variant='outlined' sx={{ mt: 3, p: { xs: 2, sm: 3 }, borderRadius: 3 }}>
                        <Stack direction={{ xs: 'column', sm: 'row' }} spacing={2} alignItems={{ sm: 'center' }}>
                            <StorageIcon color='primary' sx={{ fontSize: 36, flexShrink: 0 }} />
                            <Box sx={{ flex: 1, minWidth: 0 }}>
                                <Typography variant='h6' sx={{ fontFamily: 'monospace', wordBreak: 'break-all' }}>
                                    {info.fileName}
                                </Typography>
                                <Typography variant='body2' color='textSecondary'>
                                    SQLite database / {formatBytes(info.sizeBytes)} / {info.tables.length} tables /
                                    release v{info.version}
                                </Typography>
                            </Box>
                            <Button
                                variant='contained'
                                color='primary'
                                size='large'
                                startIcon={<DownloadIcon />}
                                href='/api/dataset/download'
                                download={info.fileName}
                                sx={{ flexShrink: 0 }}
                            >
                                Download
                            </Button>
                        </Stack>

                        <Box sx={{ mt: 2, display: 'flex', alignItems: 'center', gap: 1, minWidth: 0 }}>
                            <Typography variant='caption' color='textSecondary' sx={{ flexShrink: 0 }}>
                                SHA-256
                            </Typography>
                            <Typography
                                variant='caption'
                                sx={{ fontFamily: 'monospace', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}
                            >
                                {info.sha256}
                            </Typography>
                            <Tooltip title='Copy checksum' arrow>
                                <IconButton size='small' onClick={() => copyText(info.sha256, 'checksum')}>
                                    <ContentCopyIcon fontSize='inherit' />
                                </IconButton>
                            </Tooltip>
                        </Box>

                        <Typography variant='body2' color='textSecondary' sx={{ mt: 1 }}>
                            The dataset is licensed under{' '}
                            <Link href={DATABASE_LICENSE_URL} target='_blank' rel='noopener noreferrer'>
                                CC BY 4.0
                                <ExitIcon fontSize='inherit' sx={{ ml: 0.5, verticalAlign: 'text-bottom' }} />
                            </Link>
                            . If you use it, please{' '}
                            <Link href={PUBLICATION_URL} target='_blank' rel='noopener noreferrer'>
                                cite our publication
                                <ExitIcon fontSize='inherit' sx={{ ml: 0.5, verticalAlign: 'text-bottom' }} />
                            </Link>
                            .
                        </Typography>
                    </Card>

                    <Box sx={{ mt: 5, display: 'flex', flexWrap: 'wrap', alignItems: 'center', gap: 1 }}>
                        <Typography variant='h5' sx={{ flex: 1, minWidth: 200 }}>
                            Database schema
                        </Typography>
                        <Button variant='outlined' startIcon={<DownloadIcon />} onClick={downloadSvg}>
                            SVG
                        </Button>
                        <Button variant='outlined' startIcon={<DownloadIcon />} onClick={downloadPng}>
                            PNG
                        </Button>
                    </Box>
                    <Typography variant='body2' color='textSecondary' sx={{ mt: 1, maxWidth: 760 }}>
                        Each box is a table with its columns and row count. Lines run from a foreign key to the key it
                        references. Domains and proteins have their names in the synonym tables, and
                        the two association tables link domains to their substrates and to their proteins.
                    </Typography>

                    {/* the diagram is built from the schema with every value escaped (utils/erDiagram.js) */}
                    <Box
                        sx={{
                            mt: 2,
                            border: 1,
                            borderColor: 'divider',
                            borderRadius: 3,
                            overflowX: 'auto',
                            backgroundColor: 'background.paper',
                            '& svg': { display: 'block', maxWidth: '100%', height: 'auto', minWidth: 720 },
                        }}
                        role='img'
                        aria-label={`Entity-relationship diagram of the ${info.tables.length} database tables`}
                        dangerouslySetInnerHTML={{ __html: diagram.svg }}
                    />

                    <Box sx={{ mt: 5, display: 'flex', alignItems: 'center', gap: 1 }}>
                        <Typography variant='h5' sx={{ flex: 1 }}>
                            Example query
                        </Typography>
                        <Tooltip title='Copy query' arrow>
                            <IconButton onClick={() => copyText(EXAMPLE_QUERY, 'query')}>
                                <ContentCopyIcon fontSize='small' />
                            </IconButton>
                        </Tooltip>
                    </Box>
                    <Box
                        component='pre'
                        sx={{
                            mt: 1,
                            p: 2,
                            borderRadius: 2,
                            backgroundColor: 'surface.sunken',
                            overflowX: 'auto',
                            fontSize: '0.8rem',
                        }}
                    >
                        {`sqlite3 ${info.fileName}\n\n${EXAMPLE_QUERY}`}
                    </Box>
                </>
            )}
        </Box>
    );
};

export default Dataset;
