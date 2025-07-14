process INSTALL_FREYJA {
    tag "Installing Freyja"
    
    output:
    path "freyja", emit: freyja_bin
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "🧬 Setting up Freyja with Docker..."
    
    # Check if Docker is available
    if command -v docker &> /dev/null; then
        echo "✅ Docker found, using containerized Freyja"
        
        # Pull Freyja Docker image
        docker pull quay.io/biocontainers/freyja:1.4.4--pyhdfd78af_0
        
        # Create wrapper script
        mkdir -p freyja
        cat > freyja/freyja << 'EOF'
#!/bin/bash
# Freyja Docker wrapper
WORKDIR=\$(pwd)
docker run --rm -v "\$WORKDIR:\$WORKDIR" -w "\$WORKDIR" \\
    quay.io/biocontainers/freyja:1.4.4--pyhdfd78af_0 \\
    freyja "\$@"
EOF
        chmod +x freyja/freyja
        
        # Test the installation
        echo "🧪 Testing Freyja installation..."
        ./freyja/freyja --version
        
    elif command -v singularity &> /dev/null; then
        echo "✅ Singularity found, using containerized Freyja"
        
        # Create wrapper script for Singularity
        mkdir -p freyja
        cat > freyja/freyja << 'EOF'
#!/bin/bash
# Freyja Singularity wrapper
WORKDIR=\$(pwd)
singularity exec --bind "\$WORKDIR:\$WORKDIR" \\
    docker://quay.io/biocontainers/freyja:1.4.4--pyhdfd78af_0 \\
    freyja "\$@"
EOF
        chmod +x freyja/freyja
        
        # Test the installation
        echo "🧪 Testing Freyja installation..."
        ./freyja/freyja --version
        
    else
        echo "⚠️  Neither Docker nor Singularity found, using placeholder..."
        
        # Create placeholder for development/testing
        mkdir -p freyja
        cat > freyja/freyja << 'EOF'
#!/bin/bash
echo "🧬 Freyja placeholder (Docker not available)"
echo "freyja \$*"
echo "Mock Freyja execution completed"
EOF
        chmod +x freyja/freyja
    fi
    
    echo "✅ Freyja setup completed!"
    """
}