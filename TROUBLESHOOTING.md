# Troubleshooting

## Error when creating the conda environment

Error message:
```
Solving environment: failed
CondaValueError: Malformed version string '~': invalid character(s)
```

This error results from an outdated conda installation.
Execute `conda update conda` to update conda to the latest version.