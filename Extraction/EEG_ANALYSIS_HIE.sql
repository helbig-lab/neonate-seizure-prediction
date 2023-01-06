/*
Query to extract EEG as well as demographic and diagnostic data for neonates
*/
SELECT DISTINCT
    p.pat_mrn_id,
    hsp.pat_enc_csn_id,
    hsp.contact_date - p.birth_date AS contact_age,
    vi.end_exam_dttm - p.birth_date AS age_eeg_days,
    op.order_proc_id,
    op.order_status_c,
    sex.name AS sex,
    dg.race AS patient_race,
    dg.ethnicity,
    diag.dx_name AS diagnosis,
    edg_current_icd10.code,
    p.death_date,
    eap.proc_name,
    op.ordering_date,
    op.proc_bgn_time,
    op.proc_end_time,
    vi.end_exam_dttm,
    df.cur_value_datetime,
    df.cur_value_user_id,
    df.cur_value_source,
    df.element_id,
    con.name AS sde_concept_name,
    dv.smrtdta_elem_value

FROM clarity.order_proc AS op
INNER JOIN clarity.ord_findings AS ordf
    ON op.order_proc_id = ordf.order_proc_id
INNER JOIN patient AS p
    ON op.pat_id = p.pat_id
INNER JOIN clarity.smrtdta_elem_data AS df
    ON ordf.cv_finding_id = df.record_id_numeric
INNER JOIN clarity.smrtdta_elem_value AS dv
    ON df.hlv_id = dv.hlv_id
INNER JOIN clarity_concept AS con
    ON df.element_id = con.concept_id
INNER JOIN clarity_eap AS eap
    ON op.proc_id = eap.proc_id
INNER JOIN pat_enc_hsp AS hsp
    ON hsp.pat_enc_csn_id = op.pat_enc_csn_id
LEFT JOIN zc_sex AS sex
    ON sex.internal_id = p.sex_c
LEFT JOIN research.v_demog AS dg
    ON dg.pat_id = p.pat_id
LEFT JOIN hsp_acct_dx_list AS hdxl
    ON hdxl.hsp_account_id = hsp.hsp_account_id
LEFT JOIN clarity_edg AS diag
    ON diag.dx_id = hdxl.dx_id
LEFT JOIN clarity.edg_current_icd10 AS edg_current_icd10
    ON diag.dx_id = edg_current_icd10.dx_id
LEFT JOIN v_img_study AS vi
    ON vi.order_id = op.order_proc_id


WHERE op.order_type_c = 5004
    AND (op.order_status_c IN (2, 3, 5) OR op.order_status_c IS NULL)
    AND hsp.adt_pat_class_c IN (1)
    AND hsp.admit_conf_stat_c IN (1, 4)
    AND hsp.contact_date >= to_date('2017-01-01', 'yyyy-mm-dd')
    AND (vi.end_exam_dttm - p.birth_date) <= 30 --filter for only patients who are under 30 days old when given EEG to reduce accidental capture of patients cooled in the past
