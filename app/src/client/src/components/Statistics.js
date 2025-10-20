import React, { useEffect, useMemo, useState } from 'react';
import { Box, Card, CardContent, CardHeader, LinearProgress, Typography } from '@mui/material';
import { ResponsiveContainer, PieChart, Pie, Tooltip, Legend, Cell } from 'recharts';
import { toast } from 'react-toastify';

const SQL_ANNOTATED_BY_DOMAIN = `
SELECT
  t.domain AS label,
  COUNT(DISTINCT ad.id) AS value
FROM adenylation_domain ad
JOIN substrate_domain_association sda ON sda.domain_id = ad.id          -- annotated only
JOIN protein_domain_association   pda ON pda.domain_id = ad.id
JOIN protein                      p   ON p.id = pda.protein_id
JOIN taxonomy                     t   ON t.id = p.taxonomy_id
GROUP BY t.domain
ORDER BY value DESC
`;

const COLORS = [
  '#8884d8', '#82ca9d', '#ffc658', '#8dd1e1', '#a4de6c',
  '#d0ed57', '#ffa07a', '#90ee90', '#87cefa', '#ffb6c1'
];

export default function Statistics() {
  const [rows, setRows] = useState([]);
  const [loading, setLoading] = useState(false);

  // Fetch once on mount
  useEffect(() => {
    let cancelled = false;
    (async () => {
      setLoading(true);
      try {
        const res = await fetch('/api/sql', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            query: SQL_ANNOTATED_BY_DOMAIN,
            page: 0,
            pageSize: 1000,
            sortBy: null,
            sortDir: null,
          }),
        });
        if (!res.ok) throw new Error(await res.text());
        const data = await res.json();
        if (!cancelled) setRows(Array.isArray(data.rows) ? data.rows : []);
      } catch (e) {
        if (!cancelled) toast.error(`Failed to load stats: ${e.message}`);
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, []);

  const chartData = useMemo(() => {
    // Expecting rows like: [{ label: 'Bacteria', value: 123 }, ...]
    return rows.map((r) => ({
      label: `${r.label} (${r.value})` ?? 'Unknown',
      value: Number(r.value ?? 0),
    })).filter(d => d.value > 0);
  }, [rows]);

  return (
    <Card sx={{ height: 420 }}>
      <CardHeader
        title="Annotated A-domains per taxonomic domain"
        subheader="Counts of distinct adenylation domains that have substrate annotations"
      />
      <CardContent sx={{ height: 340, position: 'relative' }}>
        {loading && (
          <Box sx={{ position: 'absolute', left: 0, right: 0, top: 0 }}>
            <LinearProgress />
          </Box>
        )}
        {chartData.length === 0 && !loading ? (
          <Box sx={{ height: 300, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
            <Typography variant="body2" color="text.secondary">No data</Typography>
          </Box>
        ) : (
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={chartData}
                dataKey="value"
                nameKey="label"
                outerRadius="80%"
                isAnimationActive={false}
              >
                {chartData.map((entry, idx) => (
                  <Cell key={`cell-${idx}`} fill={COLORS[idx % COLORS.length]} />
                ))}
              </Pie>
              <Tooltip />
              <Legend />
            </PieChart>
          </ResponsiveContainer>
        )}
      </CardContent>
    </Card>
  );
}

