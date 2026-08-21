
# diag_manager_nml

## Overview

The `diag_manager_nml` namelist contains runtime configuration options for the diagnostic manager.

Although the namelist is defined in `diag_manager_mod`, all variables belong to the `diag_data_mod` module. See
diag_data_mod for more details.
## Variables

### `append_pelist_name`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: If true, appends the processor element list name to output filenames in the format: `<filename>.<pelist_name>.nc`.
  Useful for distinguishing output files from different processor configurations.

### `mix_snapshot_average_fields`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Controls whether snapshot (instantaneous) and time-averaged fields can coexist in the same output file. When false, snapshot and averaged fields must be separated into different files to avoid timestamp conflicts. This option is only applicable for the legacy diag manager. The modern diag manager does not allow mixing snapshot (instantaneous) and time-averaged fields.

### `max_files`
- Type: `INTEGER`
- Default: `31`
- Description: Maximum number of output files that can be managed by the diagnostic manager. This option is only applicable for the legacy diag manager.

### `max_output_fields`
- Type: `INTEGER`
- Default: `300`
- Description: Maximum number of output fields (ie. diagnostic variables as defined in the diag_table) that can be registered with the diagnostic manager. This option is only applicable for the legacy diag manager.

### `max_input_fields`
- Type: `INTEGER`
- Default: `300`
- Description: Maximum number of input fields (ie. fields registered at runtime) that can be processed by the diagnostic manager. This option is only applicable for the legacy diag manager.

### `max_axes`
- Type: `INTEGER`
- Default: `60`
- Description: Maximum number of coordinate axes that can be registered for diagnostic fields.

### `do_diag_field_log`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Enables logging of diagnostic field registration and data sending operations. Helpful for debugging diagnostic setup and data flow.
- Notes: Logfiles are written to `diag_field_log.out.<root pe number>`. The `field_log_separator` variable controls the separator character for log entries (default is `|`).

### `write_bytes_in_file`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: If true, writes the number of bytes written to each output file as metadata. Provides information about file sizes and I/O volume. This option is only applicable for the legacy diag manager.

### `debug_diag_manager`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Enables debug mode for the diagnostic manager, providing additional diagnostic output and validation checks during execution.

### `max_num_axis_sets`
- Type: `INTEGER`
- Default: `25`
- Description: Maximum number of axis sets that can be defined. This option is only applicable for the legacy diag manager.

### `use_cmor`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Overrides user-specified missing values and instead use CMOR (Climate Model Output Rewriter) standard missing values (`-1.0e20`). Required for CMIP-compliant output.

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
- Description: Controls whether the file start date is prepended to output filenames. For example, "00010101.<filename>.nc"
- Notes: Requires that `diag_manager_init` be called with the `time_init` parameter.

### `region_out_use_alt_value`
- Type: `LOGICAL`
- Default: `.TRUE.`
- Description: Determines which sentinel value to use when checking regional output boundaries. Uses `GLO_REG_VAL_ALT` (`-1`) when true, and `GLO_REG_VAL` (`-999`) when false. This option is only applicable for the legacy diag manager. The modern diag manager can only accept -999 as an option"

### `use_modern_diag`
- Type: `LOGICAL`
- Default: `.FALSE.`
- Description: Enables the modern diagnostic manager implementation with YAML-based diag tables. When false, uses the legacy ASCII-based diagnostic manager for backward compatibility.

### `use_clock_average`
- Type: `LOGICAL`
- Default: `.FALSE.`
Description: Controls how averaging windows are defined. When true, averaging of variables is done based on the clock. For example, if enabled and you start at day 1 hour 5, a 1 day freqency will only account for the rest of the hours in that day, so 19 hours total. Normally, the averaging would be done over the full 24 hour period regardless if it goes into the next day.
