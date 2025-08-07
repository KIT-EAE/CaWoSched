import os
import glob
import re

def generate_experiments_yaml():
    # Base directories
    scheduling_graphs_dir = "SchedulingGraphs"
    profiles_base_dir = "Profiles"
    mappings_dir = "Mappings"
    
    # Output YAML content
    yaml_content = []
    yaml_content.append("instdir: \"./../\"")
    yaml_content.append("instances:")
    yaml_content.append("  - repo: local")
    yaml_content.append("    items:")
    
    # Get all .dot files from SchedulingGraphs
    dot_files = glob.glob(os.path.join(scheduling_graphs_dir, "*.dot"))
    
    for dot_file in sorted(dot_files):
        filename = os.path.basename(dot_file)
        
        # Extract the prefix and pattern (e.g., 'atacseq1000' and 's1') before '.dot'
        match = re.search(r'^(.+?)(s\d+)\.dot$', filename)
        if not match:
            continue
            
        prefix = match.group(1)  # e.g., 'atacseq1000'
        pattern = match.group(2)  # e.g., 's1', 's2'
        base_name = filename.replace('.dot', '')  # e.g., 'atacseq1000s1'
        
        # Determine size and mapping based on pattern
        if pattern == 's1':
            size = "small"
            mapping_file = f"{prefix}_mapping_small.txt"
        elif pattern == 's2':
            size = "large" 
            mapping_file = f"{prefix}_mapping_large.txt"
        else:
            # Default case for other patterns
            size = "small"
            mapping_file = f"{prefix}_mapping_small.txt"
        
        # Find corresponding profile directory
        profile_dir = os.path.join(profiles_base_dir, f"{prefix}_{size}")
        
        if not os.path.exists(profile_dir):
            print(f"Warning: Profile directory {profile_dir} not found for {filename}")
            continue
        
        # Get all profile files for this pattern
        profile_pattern = os.path.join(profile_dir, f"setting3_*.input")
        profile_files = glob.glob(profile_pattern)
        
        for profile_file in sorted(profile_files):
            profile_filename = os.path.basename(profile_file)
            
            # Extract the name part from profile filename (remove setting3_ and .input)
            name_match = re.search(r'setting3_(.+)\.input$', profile_filename)
            if not name_match:
                continue
                
            instance_name = name_match.group(1)
            
            # Generate YAML entry
            yaml_content.append(f"      - name: {instance_name}")
            yaml_content.append("        files:")
            yaml_content.append(f"          - {os.path.join(scheduling_graphs_dir, filename)}")
            yaml_content.append(f"          - {os.path.join(mappings_dir, mapping_file)}")
            yaml_content.append(f"          - {profile_file}")
            yaml_content.append("")

        # Get all profile files for this pattern
        profile_pattern = os.path.join(profile_dir, f"setting1_*.input")
        profile_files = glob.glob(profile_pattern)
        
        for profile_file in sorted(profile_files):
            profile_filename = os.path.basename(profile_file)

            # Extract the name part from profile filename (remove setting1_ and .input)
            name_match = re.search(r'setting1_(.+)\.input$', profile_filename)
            if not name_match:
                continue
                
            instance_name = name_match.group(1)
            
            # Generate YAML entry
            yaml_content.append(f"      - name: {instance_name}")
            yaml_content.append("        files:")
            yaml_content.append(f"          - {os.path.join(scheduling_graphs_dir, filename)}")
            yaml_content.append(f"          - {os.path.join(mappings_dir, mapping_file)}")
            yaml_content.append(f"          - {profile_file}")
            yaml_content.append("")

        # Get all profile files for this pattern
        profile_pattern = os.path.join(profile_dir, f"setting2_*.input")
        profile_files = glob.glob(profile_pattern)
        
        for profile_file in sorted(profile_files):
            profile_filename = os.path.basename(profile_file)

            # Extract the name part from profile filename (remove setting2_ and .input)
            name_match = re.search(r'setting2_(.+)\.input$', profile_filename)
            if not name_match:
                continue
                
            instance_name = name_match.group(1)
            
            # Generate YAML entry
            yaml_content.append(f"      - name: {instance_name}")
            yaml_content.append("        files:")
            yaml_content.append(f"          - {os.path.join(scheduling_graphs_dir, filename)}")
            yaml_content.append(f"          - {os.path.join(mappings_dir, mapping_file)}")
            yaml_content.append(f"          - {profile_file}")
            yaml_content.append("")

        # Get all profile files for this pattern
        profile_pattern = os.path.join(profile_dir, f"setting4_*.input")
        profile_files = glob.glob(profile_pattern)
        
        for profile_file in sorted(profile_files):
            profile_filename = os.path.basename(profile_file)

            # Extract the name part from profile filename (remove setting4_ and .input)
            name_match = re.search(r'setting4_(.+)\.input$', profile_filename)
            if not name_match:
                continue
                
            instance_name = name_match.group(1)
            
            # Generate YAML entry
            yaml_content.append(f"      - name: {instance_name}")
            yaml_content.append("        files:")
            yaml_content.append(f"          - {os.path.join(scheduling_graphs_dir, filename)}")
            yaml_content.append(f"          - {os.path.join(mappings_dir, mapping_file)}")
            yaml_content.append(f"          - {profile_file}")
            yaml_content.append("")

    
    # Add experiments section
    yaml_content.extend([
        "",
        "",
        "experiments:",
        "  - name: 'CaWoSched-LS'",
        "    exclusive: true",
        "    stdout: out",
        "    args: ['./../../../deploy/multi_machine_scheduler', '@INSTANCE:0@', '@INSTANCE:1@', '@INSTANCE:2@']",
        "",
        "  - name: 'baseline'",
        "    exclusive: true", 
        "    stdout: out",
        "    args: ['./../../../deploy/multi_machine_scheduler', '@INSTANCE:0@', '@INSTANCE:1@', '@INSTANCE:2@', '--baseline_only']",
        "",
        "  - name: 'CaWoSched-no-LS'",
        "    exclusive: true",
        "    stdout: out", 
        "    args: ['./../../../deploy/multi_machine_scheduler', '@INSTANCE:0@', '@INSTANCE:1@', '@INSTANCE:2@', '--no_LS']"
    ])
    
    return "\n".join(yaml_content)

if __name__ == "__main__":
    # Generate the YAML content
    yaml_output = generate_experiments_yaml()
    
    # Write to experiments.yml
    with open("experiments.yml", "w") as f:
        f.write(yaml_output)
    
    print("Generated experiments.yml successfully!")
    print("\nFirst few lines of output:")
    print("\n".join(yaml_output.split("\n")[:20]) + "\n...")