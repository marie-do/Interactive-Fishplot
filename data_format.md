## 2.1 The Data

The development of the customized fishplot function is based on the Morita et al. (2020) dataset.
This dataset includes 153 acute myeloid leukemia (AML) point mutation trees from 123 patients,
obtained using scDNA-seq. Tree inference was performed using the SCITE method.

The data are organized as a nested JSON-like structure representing clonal evolution trees with
associated metadata.

### Sample naming

Each sample is named using:
AML-{PATIENT_ID}-{TIMEPOINT}

Example: AML-63-005

### Tree structure

- Each sample corresponds to a single tree
- Nodes are hierarchically linked
- Each node may have zero or more children
- The root node represents the ancestral clone

### Node attributes

| Attribute        | Type                | Description |
|------------------|---------------------|-------------|
| node_id          | string / integer    | Unique node identifier |
| label            | list (OPTIONAL)     | Typically ["Root"] |
| name             | string              | String representation |
| matching_label   | integer(OPTIONAL)   | Mapping identifier |
| size_percent     | float (OPTIONAL)    | Proportion of cells |
| gene_events      | dictionary          | Gene → mutation mapping |

### Example JSON structure for a single sample (AML-04-001):

```json
{
  "AML-04-001": {
    "tree": {
      "node_id": "12",
      "label": ["Root"],
      "variant": "",
      "name": "12",
      "matching_label": 4,
      "children": [
        {
          "node_id": "6",
          "size_percent": 0.023,
          "matching_label": 3,
          "gene_events": {
            "SF3B1": { "SNV": "p.K666N" }
          },
          "children": [
            {
              "node_id": "7",
              "size_percent": 0.04,
              "matching_label": 15,
              "gene_events": {
                "SRSF2": { "SNV": "p.P95H" }
              },
              "children": [
                {
                  "node_id": "0",
                  "size_percent": 0.701,
                  "matching_label": 6,
                  "gene_events": {
                    "FLT3-ITD": { "SNV": "" }
                  },
                  "children": [
                    {
                      "node_id": "5",
                      "size_percent": 0.175,
                      "matching_label": 17,
                      "gene_events": {
                        "PTPN11": { "SNV": "p.F71L" }
                      }
                    }
                  ]
                },
                {
                  "node_id": "1",
                  "size_percent": 0.007,
                  "matching_label": 7,
                  "gene_events": {
                    "FLT3": { "SNV": "p.D835E" }
                  },
                  "children": [
                    {
                      "node_id": "9",
                      "size_percent": 0.001,
                      "matching_label": 9,
                      "gene_events": {
                        "WT1": { "SNV": "p.R380W" }
                      },
                      "children": [
                        {
                          "node_id": "10",
                          "size_percent": 0.044,
                          "matching_label": 9,
                          "gene_events": {
                            "WT1": { "SNV": "p.R380fs" }
                          }
                        }
                      ]
                    }
                  ]
                },
                {
                  "node_id": "4",
                  "size_percent": 0.009,
                  "matching_label": 14,
                  "gene_events": {
                    "NRAS": { "SNV": "p.G12D" }
                  },
                  "children": [
                    {
                      "node_id": "2",
                      "size_percent": 0.001,
                      "matching_label": 1,
                      "gene_events": {
                        "IDH1": { "SNV": "p.R132C" }
                      }
                    }
                  ]
                }
              ]
            }
          ]
        }
      ]
    },
    "metadata": {
      "VitalStatus": "Dead NOS"
    }
  }
}
```

### To add a new patient:
Create a new entry in the JSON structure with the same format as above, ensuring that the tree structure and node attributes are consistent with the existing data.
