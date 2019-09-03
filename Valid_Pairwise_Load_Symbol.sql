CALL create_staging_tables();
-- Setup Valid Load
DELETE FROM GENE WHERE id > 0;
INSERT INTO GENE (symbol, external_id, external_id_source, description)
VALUES ("Gene1", "GN001", NULL, NULL);

INSERT INTO GENE (symbol, external_id, external_id_source, description)
VALUES ("Gene2", "GN002", NULL, NULL);

INSERT INTO GENE (symbol, external_id, external_id_source, description)
VALUES ("Gene3", "GN003", NULL, NULL);

INSERT INTO GENE (symbol, external_id, external_id_source, description)
VALUES ("Gene4", "GN004", NULL, NULL);

INSERT INTO stg_EXPERIMENT (lab, title, organism)
VALUES ("DNA", "Trial 1", "Mouse");

INSERT INTO stg_SAMPLE (label, external_accession, platform, proc_batch)
VALUES ("S1", NULL, "RNA-Seq", 1);
INSERT INTO stg_SAMPLE (label, external_accession, platform, proc_batch)
VALUES ("S2", NULL, "RNA-Seq", 1);
INSERT INTO stg_SAMPLE (label, external_accession, platform, proc_batch)
VALUES ("S3", NULL, "RNA-Seq", 1);
INSERT INTO stg_SAMPLE (label, external_accession, platform, proc_batch)
VALUES ("S4", NULL, "RNA-Seq", 1);
INSERT INTO stg_SAMPLE (label, external_accession, platform, proc_batch)
VALUES ("S5", NULL, "RNA-Seq", 1);
INSERT INTO stg_SAMPLE (label, external_accession, platform, proc_batch)
VALUES ("S6", NULL, "RNA-Seq", 1);

INSERT INTO stg_SAMPLE_NUM_ATR (sample_label, attr_name, attr_def, attr_units, attr_value)
VALUES ("S1", "Age", "Mouse Age in Years", "Years", 1);
INSERT INTO stg_SAMPLE_NUM_ATR (sample_label, attr_name, attr_def, attr_units, attr_value)
VALUES ("S2", "Age", "Mouse Age in Years", "Years", 1);
INSERT INTO stg_SAMPLE_NUM_ATR (sample_label, attr_name, attr_def, attr_units, attr_value)
VALUES ("S3", "Age", "Mouse Age in Years", "Years", 1);

INSERT INTO stg_SAMPLE_NUM_ATR (sample_label, attr_name, attr_def, attr_units, attr_value)
VALUES ("S4", "Age", "Mouse Age in Years", "Years", 0.5);
INSERT INTO stg_SAMPLE_NUM_ATR (sample_label, attr_name, attr_def, attr_units, attr_value)
VALUES ("S5", "Age", "Mouse Age in Years", "Years", 0.5);
INSERT INTO stg_SAMPLE_NUM_ATR (sample_label, attr_name, attr_def, attr_units, attr_value)
VALUES ("S6", "Age", "Mouse Age in Years", "Years", 0.5);

INSERT INTO stg_SAMPLE_CAT_ATR (sample_label, attr_name, attr_def, attr_value)
VALUES ("S1", "Genotype", "The genotype of this mouse", "WT");
INSERT INTO stg_SAMPLE_CAT_ATR (sample_label, attr_name, attr_def, attr_value)
VALUES ("S3", "Genotype", "The genotype of this mouse", "WT");
INSERT INTO stg_SAMPLE_CAT_ATR (sample_label, attr_name, attr_def, attr_value)
VALUES ("S5", "Genotype", "The genotype of this mouse", "WT");

INSERT INTO stg_SAMPLE_CAT_ATR (sample_label, attr_name, attr_def, attr_value)
VALUES ("S2", "Genotype", "The genotype of this mouse", "FN");
INSERT INTO stg_SAMPLE_CAT_ATR (sample_label, attr_name, attr_def, attr_value)
VALUES ("S4", "Genotype", "The genotype of this mouse", "FN");
INSERT INTO stg_SAMPLE_CAT_ATR (sample_label, attr_name, attr_def, attr_value)
VALUES ("S6", "Genotype", "The genotype of this mouse", "FN");

INSERT INTO stg_GROUP (label, title, criteria)
VALUES ("GR1", "Young Mice", "Age < 1");
INSERT INTO stg_GROUP (label, title, criteria)
VALUES ("GR2", "Old Mice", "Age >= 1");

INSERT INTO stg_GROUP (label, title, criteria)
VALUES ("GR3", "Wildtype", "Genotype = 'WT'");
INSERT INTO stg_GROUP (label, title, criteria)
VALUES ("GR4", "Mutant", "Genotype = 'FN'");

INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S1", "GR1");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S2", "GR1");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S3", "GR1");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S4", "GR2");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S5", "GR2");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S6", "GR2");

INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S1", "GR3");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S3", "GR3");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S5", "GR3");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S2", "GR4");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S4", "GR4");
INSERT INTO stg_GROUP_ASSIGN (sample_label, group_label)
VALUES ("S6", "GR4");

INSERT INTO stg_CONTRAST (label, group_1_label, group_2_label)
VALUES ("OLDvsYOUNG", "GR1", "GR2");

INSERT INTO stg_CONTRAST (label, group_1_label, group_2_label)
VALUES ("WTvsFN", "GR3", "GR4");

INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene1", "OLDvsYOUNG", 10, 5, -1, 0.00005, 0.005);

INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene2", "OLDvsYOUNG", 10, 2, -5, 0.00005, 0.005);

INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene3", "OLDvsYOUNG", 10, 20, 2, 0.00005, 0.005);

INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene4", "OLDvsYOUNG", 10, 40, 4, 0.00005, 0.005);


INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene1", "WTvsFN", 10, 40, 4, 0.00005, 0.005);
INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene2", "WTvsFN", 10, 20, 2, 0.00005, 0.005);
INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene3", "WTvsFN", 10, 2, -5, 0.00005, 0.005);
INSERT INTO stg_PAIRWISE_RESULT (
    gene_id, contrast_label, group_1_avg, 
    group_2_avg, logFC, p_value, fdr
)
VALUES ("Gene4", "WTvsFN", 10, 5, -1, 0.00005, 0.005);


CALL load_pairwise_tables('symbol', @lr);
