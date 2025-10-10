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
    if ! command -v docker &> /dev/null; then
        echo "❌ ERROR: Docker is required but not found!"
        echo "Please install Docker to use Freyja"
        exit 1
    fi
    
    echo "✅ Docker found, setting up Freyja container"
    
    # Pull the StaPH-B Freyja Docker image
    echo "📦 Pulling StaPH-B Freyja Docker image..."
    docker pull staphb/freyja:latest
    
    # Create wrapper script directly (not in subdirectory)
    cat > freyja << 'EOF'
#!/bin/bash
# Freyja Docker wrapper using StaPH-B image
WORKDIR=\$(pwd)
docker run --rm -v "\$WORKDIR:\$WORKDIR" -w "\$WORKDIR" \\
    staphb/freyja:latest \\
    freyja "\$@"
EOF
    
    chmod +x freyja
    
    # Test the installation
    echo "🧪 Testing Freyja installation..."
    if ./freyja --version; then
        echo "✅ Freyja setup completed successfully!"
    else
        echo "❌ ERROR: Freyja test failed!"
        exit 1
    fi
    """
}