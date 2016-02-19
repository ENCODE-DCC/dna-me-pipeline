# Bismark-ENCODE-WGBS Developer Readme

<!--
TODO: Please edit this Readme.developer.md file to include information
for developers or advanced users, for example:

==> If alignment fails "Low scratch storage space" on mem3_hdd2_x8 --multi 4 then
    use mem1_hdd2_x32 --multi 6.
    NOTE: dxapp.json::ncpus is divided in halt to get --multi for bismark alignment.
    - success on dme-align-pe mem3_hdd2_x8 ncpus=8(--multi 4) 110GB of pe fastqs
    - failed  on dme-align-pe mem3_hdd2_x8 ncpus=8(--multi 4) 126GB of pe fastqs... 

-->

## Running this app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See <a
href="https://wiki.dnanexus.com/API-Specification-v1.0.0/IO-and-Run-Specifications#Run-Specification">Run
Specification</a> in the API documentation for more information about the
available instance types.
