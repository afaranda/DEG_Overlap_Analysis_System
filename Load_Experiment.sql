/*
	File: Load_Experiment.sql
    Purpose: Define Stored procedures to execute ETL from staging tables
    Created: August 15, 2019
    Author: Adam Faranda
*/

DROP PROCEDURE IF EXISTS create_staging_tables;
DELIMITER //
CREATE PROCEDURE create_staging_tables()
BEGIN
-- Create tables for staging purposes
    DROP TABLE IF EXISTS stg_EXPERIMENT;
    CREATE TABLE stg_EXPERIMENT LIKE EXPERIMENT;
    
    DROP TABLE IF EXISTS stg_SAMPLE;
    CREATE TABLE stg_SAMPLE LIKE SAMPLE;
    ALTER TABLE stg_SAMPLE DROP PRIMARY KEY;
    ALTER TABLE stg_SAMPLE DROP experiment_id;
    
    DROP TABLE IF EXISTS stg_SAMPLE_NUM_ATR;
    CREATE TABLE stg_SAMPLE_NUM_ATR LIKE SAMPLE_NUM_ATR;
    ALTER TABLE stg_SAMPLE_NUM_ATR DROP PRIMARY KEY;
    ALTER TABLE stg_SAMPLE_NUM_ATR DROP experiment_id;
    
    DROP TABLE IF EXISTS stg_SAMPLE_CAT_ATR;
    CREATE TABLE stg_SAMPLE_CAT_ATR LIKE SAMPLE_CAT_ATR;
    ALTER TABLE stg_SAMPLE_CAT_ATR DROP PRIMARY KEY;
    ALTER TABLE stg_SAMPLE_CAT_ATR DROP experiment_id;
    
    DROP TABLE IF EXISTS stg_GROUP;
    CREATE TABLE stg_GROUP LIKE `GROUP`;
    ALTER TABLE stg_GROUP DROP PRIMARY KEY;
    ALTER TABLE stg_GROUP DROP experiment_id;

    DROP TABLE IF EXISTS stg_GROUP_ASSIGN;
    CREATE TABLE stg_GROUP_ASSIGN LIKE GROUP_ASSIGN;
    ALTER TABLE stg_GROUP_ASSIGN DROP PRIMARY KEY;
    ALTER TABLE stg_GROUP_ASSIGN DROP experiment_id;
    
    DROP TABLE IF EXISTS stg_CONTRAST;
    CREATE TABLE stg_CONTRAST LIKE CONTRAST;
    ALTER TABLE stg_CONTRAST DROP PRIMARY KEY;
    ALTER TABLE stg_CONTRAST DROP experiment_id;
    
    DROP TABLE IF EXISTS stg_PAIRWISE_RESULT;
    CREATE TABLE stg_PAIRWISE_RESULT LIKE PAIRWISE_RESULT;
    ALTER TABLE stg_PAIRWISE_RESULT DROP PRIMARY KEY;
    ALTER TABLE stg_PAIRWISE_RESULT DROP experiment_id;
    ALTER TABLE stg_PAIRWISE_RESULT MODIFY gene_id VARCHAR(500);
    
    DROP TABLE IF EXISTS stg_SAMPLE_MEASUREMENT;
    CREATE TABLE stg_SAMPLE_MEASUREMENT LIKE SAMPLE_MEASUREMENT;
    ALTER TABLE stg_SAMPLE_MEASUREMENT DROP PRIMARY KEY;
    ALTER TABLE stg_SAMPLE_MEASUREMENT DROP experiment_id;
    ALTER TABLE stg_SAMPLE_MEASUREMENT MODIFY gene_id VARCHAR(500);
    
    -- Create Views to validate primary /  foreign keys in staging tables
    -- stg_EXPERIMENT Must Only have one row
    CREATE OR REPLACE VIEW vw_primary_key_check AS
    SELECT "stg_EXPERIMENT" AS `Table`, COUNT(*) - 1 AS Check_sum  FROM stg_EXPERIMENT
        UNION ALL
    SELECT "stg_SAMPLE", COUNT(*) - COUNT(DISTINCT label)
    FROM stg_SAMPLE
        UNION ALL
    SELECT "stg_SAMPLE_CAT_ATR", COUNT(*) - COUNT(DISTINCT sample_label, attr_name)
    FROM stg_SAMPLE_CAT_ATR
        UNION ALL
    SELECT "stg_SAMPLE_NUM_ATR", COUNT(*) - COUNT(DISTINCT sample_label, attr_name)
    FROM stg_SAMPLE_NUM_ATR
        UNION ALL
    SELECT "stg_GROUP",COUNT(*) - COUNT(DISTINCT label)
    FROM stg_GROUP
        UNION ALL
    SELECT "stg_GROUP_ASSIGN", COUNT(*) - COUNT(DISTINCT sample_label, group_label) 
    FROM stg_GROUP_ASSIGN
        UNION ALL
    SELECT "stg_CONTRAST", COUNT(*) - COUNT(DISTINCT label) 
    FROM stg_CONTRAST
        UNION ALL
    SELECT "stg_PAIRWISE_RESULT", COUNT(*) - COUNT(DISTINCT contrast_label, gene_id) 
    FROM stg_PAIRWISE_RESULT
        UNION ALL
    SELECT "stg_SAMPLE_MEASUREMENT", COUNT(*) - COUNT(DISTINCT sample_label, gene_id)
    FROM stg_SAMPLE_MEASUREMENT;
    
    CREATE OR REPLACE VIEW vw_foreign_key_check AS
    SELECT 
        "stg_SAMPLE_CAT_ATR" AS `Table`,
        SUM(ISNULL(b.label)) AS Check_sum
    FROM stg_SAMPLE_CAT_ATR a
    LEFT OUTER JOIN stg_SAMPLE b
    ON a.sample_label = b.label
        UNION ALL
    SELECT
        "stg_GROUP_ASSIGN",
        SUM(ISNULL(b.label)) +
        SUM(ISNULL(c.label)) AS Check_sum
    FROM stg_GROUP_ASSIGN a
    LEFT OUTER JOIN stg_GROUP b
    ON a.group_label = b.label
    LEFT OUTER JOIN stg_SAMPLE c
    ON a.sample_label = c.label
        UNION ALL
    SELECT
        "stg_CONTRAST",
        SUM(ISNULL(b.label))
    FROM stg_CONTRAST a
    LEFT OUTER JOIN stg_GROUP b
    ON a.group_1_label = b.label
    LEFT OUTER JOIN stg_GROUP c
    ON a.group_2_label = c.label
        UNION ALL
    SELECT 
        "Invalid Contrasts", 
        COUNT(*) 
    FROM stg_GROUP_ASSIGN a
    INNER JOIN stg_CONTRAST b
    ON a.group_label = b.group_1_label
    INNER JOIN stg_GROUP_ASSIGN c
    ON c.sample_label = a.sample_label AND
    b.group_2_label = c.group_label
        UNION ALL
    SELECT
        "stg_PAIRWISE_RESULT",
        IF(COUNT(b.label) > 0, SUM(ISNULL(b.label)),0)
    FROM stg_PAIRWISE_RESULT a
    LEFT OUTER JOIN stg_CONTRAST b
    ON a.contrast_label = b.label
        UNION ALL
    SELECT
        "stg_SAMPLE_MEASUREMENT",
        IF(COUNT(b.label) > 0, SUM(ISNULL(b.label)),0)
    FROM stg_SAMPLE_MEASUREMENT a
    LEFT OUTER JOIN stg_SAMPLE b
    ON a.sample_label = b.label;
    
END //
DELIMITER ;

/* Create A procedure generating a view that reports a the list of entries
   from, stg_PAIRWISE_RESULT that do not match entries in table GENE
*/

DROP PROCEDURE IF EXISTS pairwise_nogene;
DELIMITER //
CREATE PROCEDURE pairwise_nogene(
    IN gene_id_target VARCHAR(2000)
)
pairwise:BEGIN
    
    SET @query = "CREATE OR REPLACE VIEW vw_pairwise_nogene AS";
    SET @query = CONCAT(@query, " SELECT * FROM stg_PAIRWISE_RESULT WHERE");
    SET @query = CONCAT(@query, " NOT stg_PAIRWISE_RESULT.gene_id");
    SET @query = CONCAT(@query, " IN (SELECT ", gene_id_target, " FROM GENE);");

    PREPARE stmt FROM @query;
    EXECUTE stmt;
    DEALLOCATE PREPARE stmt;
    
END //
DELIMITER ;

/* Create A procedure generating a view that reports a the list of entries
   from, stg_SAMPLE_MEASUREMENT that do not match entries in table GENE
*/

DROP PROCEDURE IF EXISTS measurement_nogene;
DELIMITER //
CREATE PROCEDURE measurement_nogene(IN gene_id_target VARCHAR(2000))
measurement:BEGIN

    SET @query = "CREATE OR REPLACE VIEW vw_measurement_nogene AS";
    SET @query = CONCAT(@query, " SELECT * FROM stg_SAMPLE_MEASUREMENT WHERE");
    SET @query = CONCAT(@query, " NOT stg_SAMPLE_MEASUREMENT.gene_id");
    SET @query = CONCAT(@query, " IN (SELECT ", gene_id_target, " FROM GENE);");

    PREPARE stmt FROM @query;
    EXECUTE stmt;
    DEALLOCATE PREPARE stmt;
    
END //
DELIMITER ;

/*
    Consider adding a procdure to determines whether data instaging tables has
    already been loaded into the system.
*/

DROP PROCEDURE IF EXISTS report_counts;
DELIMITER //
CREATE PROCEDURE report_counts()
-- Consider using a prepared statement to convert this into a view. 
report:BEGIN
    DECLARE exp_id INT(100);
    SELECT id FROM stg_EXPERIMENT INTO exp_id;
    SELECT 
        "EXPERIMENT" AS Table_Name,
        (SELECT COUNT(*) FROM stg_EXPERIMENT) AS Records_Staged,
        (SELECT COUNT(*) FROM EXPERIMENT WHERE id = exp_id) AS Records_Loaded
    UNION ALL
    SELECT 
        "SAMPLE",
        (SELECT COUNT(*) FROM stg_SAMPLE),
        (SELECT COUNT(*) FROM SAMPLE WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "SAMPLE_CAT_ATR",
        (SELECT COUNT(*) FROM stg_SAMPLE_CAT_ATR),
        (SELECT COUNT(*) FROM SAMPLE_CAT_ATR WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "SAMPLE_NUM_ATR",
        (SELECT COUNT(*) FROM stg_SAMPLE_NUM_ATR),
        (SELECT COUNT(*) FROM SAMPLE_NUM_ATR WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "GROUP",
        (SELECT COUNT(*) FROM stg_GROUP),
        (SELECT COUNT(*) FROM `GROUP` WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "GROUP_ASSIGN",
        (SELECT COUNT(*) FROM stg_GROUP_ASSIGN),
        (SELECT COUNT(*) FROM GROUP_ASSIGN WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "CONTRAST",
        (SELECT COUNT(*) FROM stg_CONTRAST),
        (SELECT COUNT(*) FROM CONTRAST WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "PAIRWISE_RESULT",
        (SELECT COUNT(*) FROM stg_PAIRWISE_RESULT),
        (SELECT COUNT(*) 
            FROM PAIRWISE_RESULT 
            WHERE experiment_id = exp_id)
    UNION ALL
    SELECT 
        "SAMPLE_MEASUREMENT",
        (SELECT COUNT(*) FROM stg_SAMPLE_MEASUREMENT),
        (SELECT COUNT(*) 
            FROM SAMPLE_MEASUREMENT 
            WHERE experiment_id = exp_id);
END //
DELIMITER ;

DROP PROCEDURE IF EXISTS load_pairwise_tables;
DELIMITER //
CREATE PROCEDURE load_pairwise_tables(IN gene_id_target VARCHAR(2000), OUT load_result VARCHAR(2000))
load_pairwise:BEGIN
    
    DECLARE count_experiment INT;
    DECLARE exp_id INT;
    
    
    START TRANSACTION;
    
    SELECT COUNT(*) FROM stg_EXPERIMENT INTO count_experiment;
    IF count_experiment = 1 THEN
        SELECT id FROM EXPERIMENT WHERE
            lab = (SELECT lab FROM stg_EXPERIMENT) AND
            title = (SELECT title FROM stg_EXPERIMENT) AND
            repository = (SELECT repository FROM stg_EXPERIMENT) AND
            external_accession = (SELECT external_accession FROM stg_EXPERIMENT) AND
            organism = (SELECT organism FROM stg_EXPERIMENT) AND
            proc_method = (SELECT proc_method FROM stg_EXPERIMENT)
        INTO exp_id;
        IF ISNULL(exp_id) THEN
            INSERT INTO EXPERIMENT (
                lab, 
                title, 
                repository, 
                external_accession,
                organism,
                proc_method
            )
            SELECT 
                lab,
                title,
                repository,
                external_accession,
                organism,
                proc_method
            FROM stg_EXPERIMENT;
            
            SELECT id FROM EXPERIMENT WHERE
                lab = (SELECT lab FROM stg_EXPERIMENT) AND
                title = (SELECT title FROM stg_EXPERIMENT) AND
                repository = (SELECT repository FROM stg_EXPERIMENT) AND
                external_accession = (SELECT external_accession FROM stg_EXPERIMENT) AND
                organism = (SELECT organism FROM stg_EXPERIMENT) AND
                proc_method = (SELECT proc_method FROM stg_EXPERIMENT)
            INTO exp_id;
            UPDATE stg_EXPERIMENT
                SET id = exp_id; 
        END IF;
    ELSE
        SET load_result = 'Failed to load due to multiple experiments';
        ROLLBACK;
        LEAVE load_pairwise;
    END IF;

    IF (SELECT SUM(Check_sum) FROM vw_primary_key_check) = 0
    THEN 
        IF (SELECT SUM(Check_sum) FROM vw_foreign_key_check) = 0
        THEN IF gene_id_target = 'symbol' 
            THEN
                -- Load pairwise contrasts, match on gene symbol
                -- DELETE FROM stg_PAIRWISE_RESULT 
                -- WHERE NOT gene_id IN (SELECT DISTINCT symbol FROM GENE);
                
                DROP TEMPORARY TABLE IF EXISTS ld_pairwise;
                CREATE TEMPORARY TABLE ld_pairwise AS
                SELECT
                    exp_id AS experiment_id,
                    a.contrast_label,
                    b.id AS gene_id,
                    a.group_1_avg,
                    a.group_2_avg,
                    a.logFC,
                    a.p_value,
                    a.fdr
                FROM stg_PAIRWISE_RESULT AS a
                INNER JOIN GENE AS b
                ON a.gene_id = b.symbol;
                
                -- Load individual measurements, match on gene symbol
                -- DELETE FROM stg_SAMPLE_MEASUREMENT
                -- WHERE NOT gene_id IN (SELECT DISTINCT symbol FROM GENE);
                
                DROP TEMPORARY TABLE IF EXISTS ld_measurement;
                CREATE TEMPORARY TABLE ld_measurement AS
                SELECT
                    exp_id AS experiment_id,
                    a.sample_label,
                    b.id AS gene_id,
                    a.measurement,
                    a.meas_type
                FROM stg_SAMPLE_MEASUREMENT AS a
                INNER JOIN GENE AS b
                ON a.gene_id = b.symbol;
                
                CALL pairwise_nogene('symbol');
                CALL measurement_nogene('symbol');
                
            ELSEIF gene_id_target = 'external_id' THEN
                -- Load pairwise contrasts, match on external identifier
                -- DELETE FROM stg_PAIRWISE_RESULT 
                -- WHERE NOT gene_id IN (SELECT DISTINCT external_id FROM GENE);
                
                DROP TEMPORARY TABLE IF EXISTS ld_pairwise;
                CREATE TEMPORARY TABLE ld_pairwise AS
                SELECT
                    exp_id AS experiment_id,
                    a.contrast_label,
                    b.id AS gene_id,
                    a.group_1_avg,
                    a.group_2_avg,
                    a.logFC,
                    a.p_value,
                    a.fdr
                FROM stg_PAIRWISE_RESULT AS a
                INNER JOIN GENE AS b
                ON a.gene_id = b.external_id
                WHERE a.gene_id IN (SELECT DISTINCT external_id FROM GENE);
                
                -- Load individual measurements, match on external_id
                -- DELETE FROM stg_SAMPLE_MEASUREMENT
                -- WHERE NOT gene_id IN (SELECT DISTINCT external_id FROM GENE);
                
                DROP TEMPORARY TABLE IF EXISTS ld_measurement;
                CREATE TEMPORARY TABLE ld_measurement AS
                SELECT
                    exp_id AS experiment_id,
                    a.sample_label,
                    b.id AS gene_id,
                    a.measurement,
                    a.meas_type
                FROM stg_SAMPLE_MEASUREMENT AS a
                INNER JOIN GENE AS b
                ON a.gene_id = b.external_id;
                
                CALL pairwise_nogene('external_id');
                CALL measurement_nogene('external_id');

            ELSE
                SET load_result = 'failed due to error in gene id target match';
                ROLLBACK;
                LEAVE load_pairwise;
            END IF;
            
            INSERT INTO SAMPLE (
                experiment_id,
                label,
                external_accession,
                short_name,
                platform,
                proc_batch
            )
            SELECT
                exp_id AS experiment_id,
                label,
                external_accession,
                short_name,
                platform,
                proc_batch
            FROM stg_SAMPLE;
            
            -- select from stg_sample_num_atr into SAMPLE_NUM_ATR
            INSERT INTO SAMPLE_NUM_ATR (
                experiment_id,
                sample_label,
                attr_name,
                attr_def,
                attr_units,
                attr_value
            )
            SELECT
                exp_id AS experiment_id,
                sample_label,
                attr_name,
                attr_def,
                attr_units,
                attr_value
            FROM stg_SAMPLE_NUM_ATR;
            
            -- select from stg_sample_cat_atr into SAMPLE_CAT_ATR
            INSERT INTO SAMPLE_CAT_ATR (
                experiment_id,
                sample_label,
                attr_name,
                attr_def,
                attr_value
            )
            SELECT
                exp_id AS experiment_id,
                sample_label,
                attr_name,
                attr_def,
                attr_value
            FROM stg_SAMPLE_CAT_ATR;
            
            -- select from stg_group into GROUP    
            INSERT INTO `GROUP` (
                experiment_id,
                label,
                title,
                criteria
            )
            SELECT
                exp_id AS experiment_id,
                label,
                title,
                criteria
            FROM stg_GROUP;
            
             -- select from stg_group_assign into GROUP_ASSIGN   
            INSERT INTO `GROUP_ASSIGN` (
                experiment_id,
                sample_label,
                group_label
            )
            SELECT
                exp_id AS experiment_id,
                sample_label,
                group_label
            FROM stg_GROUP_ASSIGN;
    
            -- select from stg_contrast into CONTRAST
            INSERT INTO CONTRAST (
                experiment_id,
                label,
                title,
                group_1_label,
                group_2_label
            )
            SELECT
                exp_id AS experiment_id,
                label,
                title,
                group_1_label,
                group_2_label
            FROM stg_CONTRAST;
        -- load PAIRWISE RESULT Table    
        INSERT INTO PAIRWISE_RESULT (
            gene_id, 
            experiment_id,
            contrast_label,
            group_1_avg,
            group_2_avg,
            logFC,
            p_value,
            fdr
        )
        SELECT
            gene_id, 
            experiment_id,
            contrast_label,
            group_1_avg,
            group_2_avg,
            logFC,
            p_value,
            fdr
        FROM ld_pairwise;

        INSERT INTO SAMPLE_MEASUREMENT (
            experiment_id,
            sample_label,
            gene_id, 
            measurement,
            meas_type
        )
        SELECT
            experiment_id,
            sample_label,
            gene_id,
            measurement,
            meas_type
        FROM ld_measurement;
            
        ELSE
            SET load_result = 'failed due to foreign key errors see vw_foreign_key_check';
            ROLLBACK;
            LEAVE load_pairwise;
        END IF;
    ELSE
        SET load_result = 'failed due to primary key errors see vw_primary_key_check';
        ROLLBACK;
        LEAVE load_pairwise;
    END IF;
    
-- Clear temporary tables, close transaction
   
            
    DROP TEMPORARY TABLE IF EXISTS ld_pairwise;
    DROP TEMPORARY TABLE IF EXISTS ld_measurement;
    SET load_result="load successful";
    COMMIT;
END //
DELIMITER ;



DROP PROCEDURE IF EXISTS drop_staging_tables;
DELIMITER //
CREATE PROCEDURE drop_staging_tables()
BEGIN
    START TRANSACTION;
    DROP VIEW IF EXISTS vw_primary_key_check;
    DROP VIEW IF EXISTS vw_foreign_key_check;
    DROP VIEW IF EXISTS vw_pairwise_nogene;
    DROP VIEW IF EXISTS vw_measurement_nogene;
    DROP TABLE IF EXISTS stg_EXPERIMENT;
    DROP TABLE IF EXISTS stg_SAMPLE;
    DROP TABLE IF EXISTS stg_SAMPLE_NUM_ATR;
    DROP TABLE IF EXISTS stg_SAMPLE_CAT_ATR;
    DROP TABLE IF EXISTS stg_GROUP;
    DROP TABLE IF EXISTS stg_GROUP_ASSIGN;
    DROP TABLE IF EXISTS stg_CONTRAST;
    DROP TABLE IF EXISTS stg_PAIRWISE_RESULT;
    DROP TABLE IF EXISTS stg_SAMPLE_MEASUREMENT;
    COMMIT;
END //
DELIMITER ;