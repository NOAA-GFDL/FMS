
# diag_manager_nml

## Overview

The `diag_manager_nml` namelist contains runtime configuration options for the diagnostic manager.

Although the namelist is defined in `diag_manager_mod`, all variables belong to the `diag_data_mod` module.
## Variables

### `append_pelist_name`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: If true, appends the processor element list name to output filenames. Useful for distinguishing output files from different processor configurations.

### `mix_snapshot_average_fields`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Controls whether snapshot (instantaneous) and time-averaged fields can coexist in the same output file. When false, snapshot and averaged fields must be separated into different files to avoid timestamp conflicts.

### `max_files`
- Type: `INTEGER`
- Default: `31`
- Description: Maximum number of output files that can be managed simultaneously by the diagnostic manager.

### `max_output_fields`
- Type: `INTEGER`
- Default: `300`
- Description: Maximum number of output fields that can be registered with the diagnostic manager.

### `max_input_fields`
- Type: `INTEGER`
- Default: `300`
- Description: Maximum number of input fields that can be processed by the diagnostic manager.

### `max_axes`
- Type: `INTEGER`
- Default: `60`
- Description: Maximum number of coordinate axes that can be defined for diagnostic fields.

### `do_diag_field_log`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Enables logging of diagnostic field registration and data sending operations. Helpful for debugging diagnostic setup and data flow.
- Notes: Logfiles are written to `diag_field_log.out.<root pe number>`. The `field_log_separator` variable controls the separator character for log entries (default is `|`).

### `write_bytes_in_files`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: If true, writes the number of bytes written to each output file. Provides information about file sizes and I/O volume.

### `debug_diag_manager`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Enables debug mode for the diagnostic manager, providing additional diagnostic output and validation checks during execution.

### `max_num_axis_sets`
- Type: `INTEGER`
- Default: `25`
- Description: Maximum number of axis sets that can be defined for organizing coordinate systems.

### `use_cmor`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Forces the use of CMOR (Climate Model Output Rewriter) standard missing values (`-1.0e20`) instead of user-specified missing values. Required for CMIP-compliant output.

### `issue_oor_warnings`
- Type: `LOGICAL`
- Default: `.TRUE.`
- Description: Controls whether warnings are issued when diagnostic field values fall outside the valid range specified during field registration.

### `oor_warnings_fatal`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Determines whether out-of-range warnings should be treated as fatal errors, causing the model to abort when invalid values are detected.

### `max_field_attributes`
- Type: `INTEGER`
- Default: `4`
- Description: Maximum number of user-defined attributes that can be attached to each diagnostic field.

### `max_file_attributes`
- Type: `INTEGER`
- Default: `2`
- Description: Maximum number of user-defined global attributes that can be attached to each output file.

### `prepend_date`
- Type: `LOGICAL`
- Default: `.TRUE.`
- Description: Controls whether the file start date is prepended to output filenames.
- Notes: Requires that `diag_manager_init` be called with the `time_init` parameter.

### `region_out_use_alt_value`
- Type: `LOGICAL`
- Default: `.TRUE.`
- Description: Determines which sentinel value to use when checking regional output boundaries. Uses `GLO_REG_VAL_ALT` (`-1`) when true, and `GLO_REG_VAL` (`-999`) when false.

### `use_mpp_io`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Selects the I/O backend: true uses `mpp_io`, false uses `fms2_io` (recommended).

### `use_modern_diag`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Enables the modern diagnostic manager implementation with YAML-based diag tables. When false, uses the legacy ASCII-based diagnostic manager for backward compatibility.

### `use_clock_average`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Controls time averaging method: true uses wall-clock time, false uses sample counting. Clock-based averaging is more accurate for variable timestep models.
