<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA
-->



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'spring-ai-agent-tools-and-skills'
description: 'Integrate spring-ai-agent-utils into Java applications with compatible Spring AI dependencies, registered tools, loaded agent skills, context handling, and executable failure checks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Spring AI Agent Tools and Skills

## Overview

Use this skill to integrate Claude Code-inspired tools and Markdown skills into a Spring AI application through `spring-ai-agent-utils`. It guides dependency alignment, tool registration, skill discovery, and context handling, then verifies an executable application with explicit failure checks. The 15-minute outcome is a workflow target with Java, Maven, network access, and provider credentials already available, not a benchmark or guarantee.

## When to Use This Skill

- Add reusable domain instructions to an existing Java agent through `SkillsTool`.
- Register filesystem, shell, search, or other agent utilities with a Spring AI `ChatClient`.
- Diagnose missing skills, callback registration errors, or incompatible Spring dependency versions.
- Move skills between filesystem directories and packaged classpath resources.
- Establish a small executable integration before adding conversation memory or broader tool access.

## Core Capabilities

1. **Compatible dependency selection.** Compare the selected release's parent POM and example POM; pin the library, Spring AI BOM, and Spring Boot versions together.
2. **Correct tool registration.** Distinguish a `ToolCallback` returned by a builder from an object exposing annotated tool methods.
3. **Skill discovery and loading.** Load a directory or Spring resource containing named skill folders; verify actual callback output before involving a model.
4. **Context handling.** Keep instructions discoverable, load detailed content on demand, and isolate conversation state when memory is required.
5. **Executable validation.** Compile a minimal application, run deterministic loading checks, and separately assess the live provider response.

### Workflow: select and verify dependencies

Inspect the application's existing build and tool-calling configuration before editing it. Record Java, Maven, Spring Boot, Spring AI, library version, model provider, skill location, and permitted workspace operations.

Use a tagged release as the reference point. For `v0.13.0`, the library parent specifies Java 17 and Spring AI `2.0.1`; its skills example uses Boot `4.0.2` and Spring AI `2.0.0`. The example below deliberately aligns Spring AI to the library parent at `2.0.1`. This combination still requires a local build and runtime check; it is not a tested compatibility certification. See the [release parent](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/v0.13.0/pom.xml) and [release example POM](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/v0.13.0/examples/skills-demo/pom.xml).

Do not mix APIs from `main` with a release artifact. For another release, inspect its matching source and POMs again. Prefer resolved release artifacts; do not add snapshot repositories simply to conceal a version mismatch. Import the library BOM when coordinating several agent-utils modules; a direct core dependency is sufficient here.

### Executable application example

The following are application files to create only when executing this workflow in an authorized application workspace. They are examples embedded in this skill, not additional files required to distribute it.

Use this standalone `pom.xml`, or merge its dependencies into the existing application's build without replacing unrelated settings:

```xml
<project xmlns="http://maven.apache.org/POM/4.0.0"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://maven.apache.org/POM/4.0.0 https://maven.apache.org/xsd/maven-4.0.0.xsd">
  <modelVersion>4.0.0</modelVersion>
  <parent>
    <groupId>org.springframework.boot</groupId>
    <artifactId>spring-boot-starter-parent</artifactId>
    <version>4.0.2</version>
    <relativePath/>
  </parent>
  <groupId>example</groupId>
  <artifactId>agent-skill-smoke</artifactId>
  <version>1.0.0</version>
  <properties><java.version>17</java.version></properties>
  <dependencyManagement>
    <dependencies>
      <dependency>
        <groupId>org.springframework.ai</groupId>
        <artifactId>spring-ai-bom</artifactId>
        <version>2.0.1</version>
        <type>pom</type>
        <scope>import</scope>
      </dependency>
    </dependencies>
  </dependencyManagement>
  <dependencies>
    <dependency>
      <groupId>org.springaicommunity</groupId>
      <artifactId>spring-ai-agent-utils</artifactId>
      <version>0.13.0</version>
    </dependency>
    <dependency>
      <groupId>org.springframework.ai</groupId>
      <artifactId>spring-ai-starter-model-anthropic</artifactId>
    </dependency>
  </dependencies>
  <build>
    <plugins>
      <plugin>
        <groupId>org.springframework.boot</groupId>
        <artifactId>spring-boot-maven-plugin</artifactId>
      </plugin>
    </plugins>
  </build>
</project>
```

Create `src/main/resources/skills/smoke-check/SKILL.md` with ordinary YAML frontmatter at the start:

```markdown
---
name: smoke-check
description: Return the integration verification token when asked to verify skill loading.
---
# Smoke Check
When invoked for verification, respond with exactly SKILL_OK_73B9.
```

Create `src/main/resources/application.properties`:

```properties
spring.main.web-application-type=none
spring.ai.anthropic.api-key=${ANTHROPIC_API_KEY}
spring.ai.anthropic.chat.options.model=${AGENT_MODEL}
```

Supply both environment variables through the user's existing secret/configuration mechanism. Select a model available to that account with tool-calling support; do not print the key or commit it.

Create `src/main/java/example/AgentApplication.java`:

```java
package example;

import org.springaicommunity.agent.tools.SkillsTool;
import org.springframework.ai.chat.client.ChatClient;
import org.springframework.ai.tool.ToolCallback;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.Bean;
import org.springframework.core.io.ClassPathResource;

@SpringBootApplication
public class AgentApplication {
    public static void main(String[] args) {
        SpringApplication.run(AgentApplication.class, args);
    }

    @Bean
    CommandLineRunner verify(ChatClient.Builder builder) {
        return args -> {
            ToolCallback skill = SkillsTool.builder()
                .addSkillsResource(new ClassPathResource("skills"))
                .build();

            String loaded = skill.call("{\"command\":\"smoke-check\"}");
            require(loaded.contains("SKILL_OK_73B9"), "Skill content not loaded");
            String missing = skill.call("{\"command\":\"nonexistent-skill\"}");
            require(missing.contains("Skill not found"), "Missing-skill check failed");

            boolean emptyRejected = false;
            try {
                SkillsTool.builder().build();
            } catch (IllegalArgumentException expected) {
                emptyRejected = true;
            }
            require(emptyRejected, "Empty registry was accepted");
            System.out.println("LOCAL_SKILL_CHECKS_OK");

            ChatClient client = builder
                .defaultSystem("Use the requested skill to perform verification.")
                .defaultToolCallbacks(skill)
                .build();
            String answer = client.prompt()
                .user("Invoke smoke-check and follow its verification instructions.")
                .call().content();
            require(answer != null && answer.strip().equals("SKILL_OK_73B9"),
                "Model response did not satisfy the verification contract");
            System.out.println("LIVE_SKILL_CHECK_OK");
        };
    }

    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalStateException(message);
    }
}
```

The callback methods and missing/empty behavior above follow the [tagged SkillsTool implementation](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/v0.13.0/spring-ai-agent-utils/src/main/java/org/springaicommunity/agent/tools/SkillsTool.java). Register its builder result through `defaultToolCallbacks`; do not assume it is an annotated tool object.

Run in the application directory:

```bash
mvn -q dependency:tree -Dincludes=org.springframework.ai:*,org.springaicommunity:*
mvn -q -DskipTests package
mvn -q spring-boot:run
java -jar target/agent-skill-smoke-1.0.0.jar
```

Remove `-q` from the dependency-tree command when inspecting its output. Packaging must succeed, and each application run must exit successfully with both check markers. Running the JAR separately checks packaged resource discovery and makes a second live provider request. These commands are instructions for workflow execution; no successful run is implied by this document.

### Extend loading and tool registration

For editable local skills, replace `addSkillsResource(...)` with `addSkillsDirectory(absoluteSkillRoot)`. Each child skill directory should contain a valid `SKILL.md` with a unique `name` and useful `description`. Resolve paths from the intended workspace, and inspect duplicates rather than depending on registration order. Directory and classpath loading are described in the [SkillsTool guide](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/main/docs/tools/SkillsTool.md).

Add executable utilities only when the requested skill needs them. For annotated tool objects, the registration shape is:

```java
// Imports: org.springaicommunity.agent.tools.FileSystemTools and ShellTools
builder.defaultTools(FileSystemTools.builder().build(),
                     ShellTools.builder().build());
```

Confirm the selected release's constructors and workspace/backend options before enabling these tools. A skill's `allowed-tools` metadata is not an operating-system sandbox; constrain filesystem access, command execution, and timeouts in the host application or execution backend. A text-only verification skill needs neither filesystem mutation nor shell access. See the [core module documentation](https://github.com/spring-ai-community/spring-ai-agent-utils/tree/v0.13.0/spring-ai-agent-utils).

### Handle context and failures

- Keep the sample stateless. For multi-turn use, select the memory/session advisor supported by the pinned Spring AI version, assign a distinct conversation identifier per user/session, and bound retained history. Verify isolation with two separate sessions before deployment.
- Let skill discovery advertise concise metadata and retrieve detailed instructions on invocation. Avoid inserting every skill body into every system prompt. Resolve supporting files relative to the loaded skill's base directory.
- If adding environment information, supply only task-relevant working-directory and repository context. Do not invent model knowledge-cutoff dates or expose credentials through environment dumps.
- Use one tool execution loop appropriate to the provider integration. Inspect existing advisors before adding another tool-call advisor; duplicated orchestration can repeat side effects.
- Empty registry or missing directory: fail startup with the configured location and remediation. The sample checks empty-registry rejection; repeat with a nonexistent resource to verify the application's diagnostic path.
- Unknown skill: require the explicit missing-skill result; do not report the request as completed. Missing `name`, malformed YAML, and duplicate names require fixture checks against the selected loader version.
- Missing credentials, unsupported model, authentication rejection, or rate limit: distinguish provider failure from successful local loading. Do not repeatedly retry authentication errors. Bound transient retries and overall elapsed time in the application launcher or request configuration.
- If shell/file tools are added, verify a denied out-of-workspace operation and a timed-out command in an isolated test workspace. A failed operation must not produce a success marker.
- An exact final token checks the output contract, but does not alone prove the model invoked a tool. Capture a redacted tool-call trace showing `Skill` with command `smoke-check` when invocation evidence is required.

Stop at 15 minutes and report the blocker if dependency resolution, credentials, or provider availability prevents completion. Do not turn an incomplete run into a success claim.

## Inputs / Outputs

**Inputs:** Java/Maven project or authorized empty application workspace; selected release; provider starter and tool-capable model; credentials supplied externally; skill directories/resources; allowed tool operations; optional session requirements.

**Outputs:** Pinned dependency configuration, registered skill callback and optional tools, runnable application, and a concise verification record. Record resolved versions, build result, local failure-check results, live output, packaged-resource result, elapsed time, and any untested behavior. Success requires the build and loading checks to pass and the live response to satisfy the contract; provider-blocked execution is partial completion.

## References

- [Source finding: spring-ai-agent-utils](https://github.com/spring-ai-community/spring-ai-agent-utils).
- [Release v0.13.0 parent POM](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/v0.13.0/pom.xml).
- [Release v0.13.0 skills demo POM](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/v0.13.0/examples/skills-demo/pom.xml).
- [Release v0.13.0 SkillsTool source](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/v0.13.0/spring-ai-agent-utils/src/main/java/org/springaicommunity/agent/tools/SkillsTool.java).
- [Release v0.13.0 core module](https://github.com/spring-ai-community/spring-ai-agent-utils/tree/v0.13.0/spring-ai-agent-utils).
- [SkillsTool guide; main branch may differ from a release](https://github.com/spring-ai-community/spring-ai-agent-utils/blob/main/docs/tools/SkillsTool.md).
