/**
 * Tables behind the network page's downloads: one row per neighbour of a
 * search, and one row per placed signature.
 *
 * Cluster ids are only meaningful at the threshold they were computed at, so
 * every row carries that threshold with it.
 */

import { SIGNATURE_LENGTH } from '../components/DistanceBadge';

const columns = (...fields) => fields.map((field) => ({ field }));

export const NEIGHBOR_COLUMNS = columns(
    'query_name', 'query_signature', 'threshold', 'rank', 'neighbor_id', 'neighbor_name', 'distance',
    'identity', 'substrates', 'genus', 'kingdom', 'cluster_id', 'neighbor_extended_signature',
);

export const QUERY_COLUMNS = columns(
    'query_name', 'query_signature', 'threshold', 'nearest_distance', 'within_threshold',
    'clusters_joined', 'cluster_sizes', 'bridges_clusters', 'nearest_neighbor', 'nearest_neighbor_substrates',
);

/**
 * One row per neighbour of a single search.
 *
 * @param {object} search - the search.
 * @param {string} search.name - query name.
 * @param {string} search.signature - query extended signature.
 * @param {number} search.threshold - threshold the cluster ids belong to.
 * @param {Array<object>} search.neighbors - neighbours as the server ranks them.
 * @returns {Array<object>} - rows keyed by NEIGHBOR_COLUMNS.
 */
export function neighborRows({ name, signature, threshold, neighbors }) {
    return (neighbors || []).map((n, i) => ({
        query_name: name,
        query_signature: signature,
        threshold,
        rank: n.rank ?? i + 1,
        neighbor_id: n.id,
        neighbor_name: n.name,
        distance: n.distance,
        identity: ((SIGNATURE_LENGTH - n.distance) / SIGNATURE_LENGTH).toFixed(3),
        substrates: n.substrates,
        genus: n.genus,
        kingdom: n.kingdom,
        cluster_id: n.cluster_id,
        neighbor_extended_signature: n.extended_signature,
    }));
}

/**
 * One row per placed signature.
 *
 * @param {Array<{key: string, name: string, signature: string}>} queries - placed signatures, in order.
 * @param {Object<string, object>} placements - server placement per query key.
 * @param {number} threshold - threshold they were placed at.
 * @returns {Array<object>} - rows keyed by QUERY_COLUMNS.
 */
export function queryRows(queries, placements, threshold) {
    return queries.filter((q) => placements[q.key]).map((q) => {
        const placement = placements[q.key];
        const clusters = placement.clusters || [];
        const nearest = (placement.neighbors || [])[0];
        return {
            query_name: q.name,
            query_signature: q.signature,
            threshold,
            nearest_distance: placement.nearest_distance,
            within_threshold: placement.within_threshold,
            clusters_joined: clusters.map((c) => c.cluster_id),
            cluster_sizes: clusters.map((c) => c.size),
            bridges_clusters: clusters.length > 1 ? 'yes' : 'no',
            nearest_neighbor: nearest?.name ?? '',
            nearest_neighbor_substrates: nearest?.substrates ?? [],
        };
    });
}

/** A timestamp for export filenames, e.g. 2026-09-21-14-03-09. */
export const exportStamp = () => new Date().toISOString().slice(0, 19).replace(/[:T]/g, '-');

/** A query name made safe for a filename. */
export const fileSafe = (name) => String(name || 'query').replace(/[^A-Za-z0-9._-]+/g, '_').slice(0, 60);
