SELECT organism, AVG(gc_content) AS avg_gc
FROM sequences
GROUP BY organism;

SELECT organism, sequence_id, gc_content
FROM sequences
ORDER BY gc_content DESC;

SELECT organism, AVG(gc_content) AS avg_gc
FROM sequences
GROUP BY organism
HAVING AVG(gc_content) > 42;